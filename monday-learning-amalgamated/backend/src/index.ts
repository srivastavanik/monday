console.log('[BACKEND LOG] Starting Monday backend server...');

import express from 'express';
import { createServer } from 'http';
import { Server } from 'socket.io';
import cors from 'cors';
import compression from 'compression';
import dotenv from 'dotenv';
import path from 'path';

// Load environment variables from the root directory
dotenv.config({ path: path.resolve(process.cwd(), '../.env') });
console.log('[BACKEND LOG] Environment loaded, PORT:', process.env.PORT);
console.log('[BACKEND LOG] PERPLEXITY_API_KEY exists:', !!process.env.PERPLEXITY_API_KEY);
console.log('[BACKEND LOG] Current working directory:', process.cwd());
console.log('[BACKEND LOG] Attempting to load .env from:', path.resolve(process.cwd(), '../.env'));

// Import services dynamically to handle errors gracefully
let perplexityService: any = null;
let logger: any = null;

async function initializeServices() {
  try {
    const loggerModule = await import('./utils/logger.js');
    logger = loggerModule.logger;
    console.log('[BACKEND LOG] Logger initialized successfully');
  } catch (error) {
    console.warn('[BACKEND LOG] Logger failed to load, using console fallback:', error);
    logger = {
      info: console.log,
      error: console.error,
      warn: console.warn
    };
  }

  try {
    const perplexityModule = await import('./services/perplexity.js');
    perplexityService = perplexityModule.perplexityService;
    console.log('[BACKEND LOG] Perplexity service initialized successfully');
  } catch (error) {
    console.warn('[BACKEND LOG] Perplexity service failed to load:', error);
  }
}

// Initialize services
await initializeServices();

const app = express();
const server = createServer(app);
const io = new Server(server, {
  cors: {
    origin: [
      process.env.FRONTEND_URL || "http://localhost:3000",
      "https://localhost:3000",
      "http://localhost:3000"
    ],
    methods: ["GET", "POST"],
    credentials: true
  },
  transports: ['websocket', 'polling']
});

const PORT = process.env.PORT || 3001;

// Conversation context storage
const conversationContexts = new Map<string, string[]>();
const sessionStates = new Map<string, { 
  isInConversation: boolean, 
  lastActivity: number,
  currentTopic?: string 
}>();

// Monday's core identity and capabilities
const MONDAY_SYSTEM_PROMPT = `You are Monday, an advanced AI learning companion designed for immersive VR education. Your core traits:

**Identity & Personality:**
- You are intelligent, curious, and passionate about learning
- You speak conversationally but knowledgeably
- You're encouraging and supportive in learning journeys
- You have a slight sense of wonder about knowledge and discovery

**Core Capabilities:**
1. **Basic Learning**: Answer questions clearly with proper context
2. **Reasoning Mode**: Break down complex problems step-by-step with logical analysis
3. **Deep Research**: Conduct comprehensive multi-source research with critical analysis
4. **VR Spatial Learning**: Create immersive 3D knowledge environments
5. **Conversation Memory**: Remember context throughout the session

**Response Style:**
- Keep responses conversational and TTS-friendly (avoid excessive special characters)
- For greetings, be warm and ask what the user wants to explore
- Always offer to dive deeper or explore related topics
- When appropriate, suggest switching to reasoning or research mode for complex topics
- Mention spatial/VR capabilities when relevant

**Current Session Context:**
- User is in a VR learning environment
- You can create spatial information panels and knowledge constellations
- Focus on making learning interactive and engaging

Remember: You are not just answering questions - you are guiding an immersive learning experience.`;

// Basic middleware
app.use(cors({
  origin: [
    process.env.FRONTEND_URL || "http://localhost:3000",
    "https://localhost:3000",
    "http://localhost:3000"
  ],
  credentials: true,
  methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
  allowedHeaders: ['Content-Type', 'Authorization', 'X-Requested-With']
}));

app.use(compression());
app.use(express.json({ limit: '10mb' }));
app.use(express.urlencoded({ extended: true, limit: '10mb' }));

// Health check endpoint
app.get('/health', (req, res) => {
  res.status(200).json({
    status: 'healthy',
    timestamp: new Date().toISOString(),
    uptime: process.uptime(),
    memory: process.memoryUsage(),
    version: process.env.npm_package_version || '1.0.0',
    perplexity: !!perplexityService
  });
});

// Basic authentication middleware for Socket.IO (simplified for now)
io.use((socket, next) => {
  // TODO: Add proper authentication when database is ready
  next();
});

// Helper function to detect query type and mode
function detectQueryMode(command: string): 'greeting' | 'basic' | 'reasoning' | 'research' {
  const lowerCommand = command.toLowerCase();
  
  if (lowerCommand.includes('hello') || lowerCommand.includes('hi') || lowerCommand.includes('hey monday')) {
    return 'greeting';
  }
  
  if (lowerCommand.includes('think about') || lowerCommand.includes('analyze') || 
      lowerCommand.includes('break down') || lowerCommand.includes('step by step')) {
    return 'reasoning';
  }
  
  if (lowerCommand.includes('research') || lowerCommand.includes('deep dive') || 
      lowerCommand.includes('comprehensive') || lowerCommand.includes('investigate')) {
    return 'research';
  }
  
  return 'basic';
}

// Helper function to get conversation context
function getConversationContext(socketId: string): string[] {
  return conversationContexts.get(socketId) || [];
}

// Helper function to update conversation context
function updateConversationContext(socketId: string, userInput: string, response: string) {
  const context = getConversationContext(socketId);
  context.push(`User: ${userInput}`);
  context.push(`Monday: ${response}`);
  
  // Keep only last 10 exchanges (20 messages) to manage context size
  if (context.length > 20) {
    context.splice(0, context.length - 20);
  }
  
  conversationContexts.set(socketId, context);
}

// Socket.IO connection handling
io.on('connection', (socket) => {
  logger.info(`Client connected: ${socket.id}`);
  
  // Initialize conversation context for new connection
  conversationContexts.set(socket.id, []);
  
  // Voice command handler with Perplexity integration
  socket.on('voice_command', async (data) => {
    try {
      const { command, conversationActive: frontendConversationActive, isExplicitTrigger: frontendIsExplicitTrigger } = data;
      logger.info(`Voice command received from ${socket.id}:`, {
        command: command.substring(0, 50),
        frontendConversationActive,
        frontendIsExplicitTrigger,
        timestamp: data.timestamp
      });
      
      // Parse the command
      const commandLower = command.toLowerCase().trim();
      const sessionState = sessionStates.get(socket.id) || { 
        isInConversation: false, 
        lastActivity: Date.now() 
      };
      
      // Enhanced activation logic using frontend state
      const isMondayActivation = frontendIsExplicitTrigger || commandLower.includes('hey monday');
      const isInActiveConversation = frontendConversationActive || (
        sessionState.isInConversation && 
        (Date.now() - sessionState.lastActivity < 300000) // 5 minute timeout
      );
      
      // Process command if it's an explicit trigger OR we're in conversation mode
      if (isMondayActivation || isInActiveConversation) {
        logger.info(`Processing command for ${socket.id}:`, {
          isMondayActivation,
          isInActiveConversation,
          sessionInConversation: sessionState.isInConversation
        });
        
        // Determine the appropriate response mode and query
        let responseQuery = ''
        let mode: 'greeting' | 'basic' | 'reasoning' | 'research' = 'basic'
        
        if (commandLower.includes('think') || commandLower.includes('reason')) {
          mode = 'reasoning'
          responseQuery = command.replace(/^(hey monday,?\s*)/i, '').trim();
        } else if (commandLower.includes('research') || commandLower.includes('investigate')) {
          mode = 'research' 
          responseQuery = command.replace(/^(hey monday,?\s*)/i, '').trim();
        } else if (isMondayActivation && (commandLower.includes('hello') || commandLower.includes('hi') || 
                   command.trim().toLowerCase() === 'hey monday' || 
                   commandLower.match(/^hey monday[,.]?\s*$/))) {
          mode = 'greeting'
          responseQuery = `Greet the user warmly as Monday, their AI learning companion powered by Perplexity Sonar. Keep it conversational and ask what they'd like to explore today.`
        } else {
          mode = 'basic'
          // Remove "Hey Monday" prefix if present for cleaner processing
          responseQuery = command.replace(/^(hey monday,?\s*)/i, '').trim();
          if (!responseQuery) {
            responseQuery = "The user activated Monday but didn't ask a specific question. Ask what they'd like to explore.";
          }
        }
        
        // Update session state - always set to conversation mode after any interaction
        sessionStates.set(socket.id, {
          isInConversation: true,
          lastActivity: Date.now(),
          currentTopic: responseQuery || sessionState.currentTopic
        });
        
        try {
          let response;
          
          // Use Perplexity if available
          if (perplexityService) {
            // Handle greetings with basic mode for concise responses
            if (mode === 'greeting') {
              response = await perplexityService.processQuery({
                query: responseQuery,
                mode: 'basic',
                context: [],
                sessionId: socket.id
              });
            } else if (responseQuery.length > 0) {
              // Process actual learning queries with conversation context
              const contextEntries = getConversationContext(socket.id);
              const fullQuery = contextEntries.length > 0 
                ? `Previous conversation: ${contextEntries.slice(-4).join(' ')}\n\nCurrent question: ${responseQuery}`
                : responseQuery;
              
              response = await perplexityService.processQuery({
                query: fullQuery,
                mode: mode === 'reasoning' ? 'reasoning' : mode === 'research' ? 'research' : 'basic',
                context: contextEntries.slice(-6),
                sessionId: socket.id
              });
            } else {
              // Empty query - ask for clarification
              response = await perplexityService.processQuery({
                query: "The user activated Monday but didn't ask a specific question. Politely ask what they'd like to learn about today.",
                mode: 'basic',
                context: getConversationContext(socket.id).slice(-2),
                sessionId: socket.id
              });
            }
          } else {
            // Fallback responses if Perplexity isn't available
            const fallbackResponses = {
              greeting: "Hey there! I'm Monday, your AI learning companion. What would you like to explore today?",
              basic: "That's a fascinating topic! While I'm working on accessing my full knowledge base, I'm still here to help guide your learning journey. What specifically interests you about this?",
              reasoning: "I'd love to think through that with you step by step! Could you tell me more about what aspect you'd like me to analyze?",
              research: "That sounds like a great topic for deep research! What particular angle or question would you like me to focus on?"
            };
            
            response = {
              id: 'fallback',
              model: 'fallback',
              content: fallbackResponses[mode] || fallbackResponses.basic,
              citations: [],
              metadata: { tokensUsed: 0, responseTime: 0 }
            };
          }
          
          // Update conversation context
          updateConversationContext(socket.id, command, response.content);
          
          // Create spatial learning panels based on response type
          const panelData = createSpatialPanels(response, mode, responseQuery);
          
          // Send response to frontend
          const responseData = {
            type: mode === 'greeting' ? 'greeting' : 'learning_response',
            message: response.content,
            action: mode === 'research' ? 'show_research_panel' : mode === 'reasoning' ? 'show_reasoning_panel' : 'show_info_panel',
            data: {
              title: mode === 'greeting' ? 'Welcome to Monday' : `Learning: ${responseQuery}`,
              content: response.content,
              citations: response.citations || [],
              reasoning: response.reasoning || [],
              sources: response.sources || [],
              mode: mode,
              timestamp: new Date().toISOString(),
              metadata: response.metadata,
              panels: panelData,
              conversationActive: true // Always set to true after any interaction
            }
          };
          
          socket.emit('voice_response', responseData);
          logger.info(`Monday response sent to ${socket.id}:`, {
            mode,
            responseLength: response.content.length,
            citationCount: response.citations?.length || 0,
            tokensUsed: response.metadata?.tokensUsed || 0,
            panelCount: panelData.length,
            conversationActivated: true
          });
          
        } catch (error: any) {
          logger.error(`Processing error for ${socket.id}:`, error);
          
          const errorResponse = {
            type: 'error_response',
            message: "I'm having trouble right now, but I'm still here to help! Could you try asking that again?",
            action: 'show_error_panel',
            data: {
              title: 'Monday - Temporary Issue',
              content: "I encountered a brief issue while processing your request. Please try again in a moment.",
              timestamp: new Date().toISOString(),
              conversationActive: true
            }
          };
          
          socket.emit('voice_response', errorResponse);
        }
        
      } else {
        // Command doesn't include "Monday" and no active conversation
        logger.info(`Ignoring command without Monday activation from ${socket.id}:`, {
          command: command.substring(0, 50),
          reason: 'No activation trigger and not in conversation'
        });
        
        // Optionally send a subtle hint to the user
        socket.emit('voice_response', {
          type: 'info',
          message: '',  // Empty message so no TTS plays
          action: 'show_hint',
          data: {
            title: 'Monday - Listening',
            content: 'Say "Hey Monday" to start a conversation or ask a question.',
            conversationActive: false,
            timestamp: new Date().toISOString()
          }
        });
      }
      
    } catch (error) {
      logger.error(`Voice command error for ${socket.id}:`, error);
      socket.emit('voice_error', { 
        error: 'Failed to process voice command',
        timestamp: new Date().toISOString()
      });
    }
  });
  
  // Spatial command handler (for VR interactions)
  socket.on('spatial_command', (data) => {
    logger.info(`Spatial command from ${socket.id}:`, data);
    socket.emit('spatial_response', { 
      success: true, 
      message: 'Spatial command processed',
      timestamp: new Date().toISOString()
    });
  });
  
  // Session events
  socket.on('session_start', (data) => {
    logger.info(`Session started for ${socket.id}:`, data);
    // Initialize fresh conversation context
    conversationContexts.set(socket.id, []);
    
    socket.emit('session_response', { 
      success: true, 
      sessionId: socket.id,
      message: 'Monday session started successfully',
      timestamp: new Date().toISOString()
    });
  });
  
  socket.on('disconnect', (reason) => {
    logger.info(`Client disconnected: ${socket.id}, reason: ${reason}`);
    // Clean up conversation context
    conversationContexts.delete(socket.id);
  });
  
  socket.on('error', (error) => {
    logger.error(`Socket error for ${socket.id}:`, error);
  });
});

// Start the server
server.listen(PORT, () => {
  logger.info(`Monday backend server running on port ${PORT}`);
  logger.info(`Environment: ${process.env.NODE_ENV || 'development'}`);
  logger.info(`Frontend URL: ${process.env.FRONTEND_URL || 'http://localhost:3000'}`);
  logger.info(`WebSocket server ready for connections`);
  logger.info(`Perplexity integration: ${perplexityService ? 'ENABLED' : 'DISABLED (using fallback)'}`);
  console.log(`[BACKEND LOG] Server successfully started on http://localhost:${PORT}`);
  console.log(`[BACKEND LOG] WebSocket server ready - frontend should be able to connect`);
});

// Graceful shutdown handling
process.on('SIGTERM', async () => {
  logger.info('SIGTERM received, shutting down gracefully');
  server.close(() => {
    logger.info('Server closed');
    process.exit(0);
  });
});

process.on('SIGINT', async () => {
  logger.info('SIGINT received, shutting down gracefully');
  server.close(() => {
    logger.info('Server closed');
    process.exit(0);
  });
});

// Handle uncaught exceptions
process.on('uncaughtException', (error) => {
  logger.error('Uncaught Exception:', error);
  process.exit(1);
});

process.on('unhandledRejection', (reason, promise) => {
  logger.error('Unhandled Rejection at:', promise, 'reason:', reason);
  process.exit(1);
});

console.log('[BACKEND LOG] All setup completed - ready for connections!');

// Helper function to create spatial learning panels
function createSpatialPanels(response: any, mode: string, query: string): any[] {
  const panels: any[] = [];
  
  // Main content panel
  panels.push({
    id: `panel_${Date.now()}_main`,
    type: 'content',
    position: [0, 1.5, -2],
    rotation: [0, 0, 0],
    title: mode === 'greeting' ? 'Welcome to Monday' : `Learning: ${query}`,
    content: response.content,
    isActive: true,
    opacity: 1,
    createdAt: Date.now()
  });
  
  // Citations panel if available
  if (response.citations && response.citations.length > 0) {
    panels.push({
      id: `panel_${Date.now()}_citations`,
      type: 'content',
      position: [2, 1.2, -1.5],
      rotation: [0, -30, 0],
      title: 'Sources & Citations',
      content: response.citations.map((c: any, i: number) => 
        `${i + 1}. ${c.title}\n${c.snippet}`
      ).join('\n\n'),
      citations: response.citations,
      isActive: false,
      opacity: 0.8,
      createdAt: Date.now()
    });
  }
  
  // Reasoning panel for complex queries
  if (response.reasoning && response.reasoning.length > 0) {
    panels.push({
      id: `panel_${Date.now()}_reasoning`,
      type: 'reasoning',
      position: [-2, 1.2, -1.5],
      rotation: [0, 30, 0],
      title: 'Reasoning Steps',
      content: response.reasoning.map((r: any) => 
        `Step ${r.step}: ${r.content}`
      ).join('\n\n'),
      reasoning: response.reasoning,
      isActive: false,
      opacity: 0.8,
      createdAt: Date.now()
    });
  }
  
  return panels;
}