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

// Conversation context storage - Updated structure
interface ConversationEntry {
  role: 'user' | 'assistant';
  content: string;
  ttsContent?: string; // Optional TTS version
  timestamp: number;
}

const conversationContexts = new Map<string, ConversationEntry[]>();
const sessionStates = new Map<string, { 
  isInConversation: boolean, 
  lastActivity: number,
  currentTopic?: string,
  isFirstInteraction: boolean  // Add this to track first interaction
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

// Helper function to determine which model to use based on command
function determineQueryMode(command: string): { mode: string; model: string } {
  const lowerCommand = command.toLowerCase().trim()
  
  // Remove "hey monday" prefix for analysis
  const cleanCommand = lowerCommand.replace(/^hey monday[,\s]*/i, '').trim()
  
  if (cleanCommand.startsWith('think about') || cleanCommand.includes('think about')) {
    return { mode: 'reasoning', model: 'sonar-reasoning-pro' }
  } else if (cleanCommand.startsWith('research') || cleanCommand.includes('research')) {
    return { mode: 'research', model: 'sonar-deep-research' }
  } else if (cleanCommand.startsWith('search the web') || cleanCommand.includes('search the web')) {
    return { mode: 'web_search', model: 'sonar-pro' }
  } else {
    // Basic conversations and queries use the regular sonar model
    return { mode: 'basic', model: 'sonar' }
  }
}

// Helper function to get conversation context
function getConversationContext(socketId: string): ConversationEntry[] {
  return conversationContexts.get(socketId) || [];
}

// Helper function to update conversation context
function updateConversationContext(
  socketId: string, 
  userInput: string, 
  fullResponse: string,
  ttsResponse?: string
) {
  const context = getConversationContext(socketId);
  
  // Store full content for API, TTS version for reference
  context.push({
    role: 'user',
    content: userInput,
    timestamp: Date.now()
  });
  
  context.push({
    role: 'assistant',
    content: fullResponse, // FULL response for API context
    ttsContent: ttsResponse, // Shortened version for TTS
    timestamp: Date.now()
  });
  
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
  sessionStates.set(socket.id, {
    isInConversation: false,
    lastActivity: Date.now(),
    isFirstInteraction: true  // Initialize as true
  });
  
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
        lastActivity: Date.now(),
        isFirstInteraction: true
      };
      
      // Check for reset conversation command
      if (commandLower.includes('reset conversation') || commandLower.includes('start over')) {
        conversationContexts.delete(socket.id);
        
        socket.emit('voice_response', {
          message: "I've reset our conversation. Let's start fresh! What would you like to explore?",
          data: {
            panels: [{
              id: `panel_${Date.now()}_reset`,
              type: 'content',
              position: [0, 1.5, -2],
              rotation: [0, 0, 0],
              title: 'Conversation Reset',
              content: "I've reset our conversation. Let's start fresh! What would you like to explore?",
              isActive: true,
              opacity: 1,
              createdAt: Date.now(),
              model: 'system'
            }],
            mode: 'basic',
            model: 'system',
            query: 'reset',
            citations: [],
            reasoning: [],
            metadata: { tokensUsed: 0, responseTime: 0 }
          }
        });
        
        return;
      }
      
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
          currentTopic: responseQuery || sessionState.currentTopic,
          isFirstInteraction: false  // Set to false after first interaction
        });
        
        try {
          console.log(`[info]: Processing command for ${socket.id}:`, {
            service: 'monday-backend',
            isMondayActivation: isMondayActivation,
            isInActiveConversation: isInActiveConversation,
            sessionInConversation: sessionState.isInConversation
          });

          // Determine query mode and model
          const { mode: determinedMode, model } = determineQueryMode(command)
          console.log(`[info]: Using ${model} model for ${determinedMode} query`)

          // Prepare the query for the AI service
          let responseQuery = command;
          if (isMondayActivation) {
            responseQuery = command.replace(/^hey monday[,\s]*/i, '').trim();
          }

          // Ensure we have meaningful content for the API
          if (!responseQuery || responseQuery.length === 0) {
            if (isMondayActivation && (commandLower.includes('hello') || commandLower.includes('hi') || 
                       command.trim().toLowerCase() === 'hey monday' || 
                       commandLower.match(/^hey monday[,.]?\s*$/))) {
              responseQuery = "The user just greeted me with 'Hey Monday'. Respond warmly as Monday, their AI learning companion, and ask what they'd like to explore or learn about today.";
            } else {
              responseQuery = "The user activated Monday but didn't provide a specific question. Ask them what they'd like to explore or learn about.";
            }
          }

          // Get AI response using the determined model
          let response;
          try {
            if (determinedMode === 'reasoning') {
              // Create progress callback for streaming updates
              const progressCallback = (update: any) => {
                console.log('🔄 Reasoning progress update:', update);
                socket.emit('reasoning_progress', {
                  type: 'progress_update',
                  update: update,
                  query: responseQuery,
                  model: model,
                  timestamp: Date.now()
                });
              };
              
              response = await perplexityService.reasoningQuery(responseQuery, getConversationContext(socket.id), progressCallback);
            } else if (determinedMode === 'research') {
              // Create progress callback for streaming updates
              const progressCallback = (update: any) => {
                console.log('🔍 Research progress update:', update);
                socket.emit('research_progress', {
                  type: 'progress_update',
                  update: update,
                  query: responseQuery,
                  model: model,
                  timestamp: Date.now()
                });
              };
              
              response = await perplexityService.deepResearch(responseQuery, getConversationContext(socket.id), progressCallback);
              
              // Extract thinking process from fullContent if not already present
              if (!response.thinkingProcess && response.fullContent) {
                const thinkMatch = response.fullContent.match(/<think>([\s\S]*?)<\/think>/i);
                if (thinkMatch) {
                  response.thinkingProcess = thinkMatch[1].trim();
                  console.log('🧠 Extracted thinking process for panel:', response.thinkingProcess.substring(0, 100) + '...');
                  
                  // Clean up main panel content by removing thinking process
                  response.fullContent = response.fullContent.replace(/<think>[\s\S]*?<\/think>/i, '').trim();
                  console.log('🧹 Cleaned main panel content');
                }
              }
            } else {
              response = await perplexityService.basicQuery(responseQuery, {
                service: 'monday-backend',
                isMondayActivation: isMondayActivation,
                isInActiveConversation: isInActiveConversation,
                sessionInConversation: sessionState.isInConversation
              });
            }
          } catch (apiError: any) {
            console.warn('Perplexity API failed, using fallback response:', apiError.message);
            
            // Create fallback response
            const fallbackResponses = {
              greeting: "Hi there! I'm Monday, your AI learning companion. I'm having some connectivity issues right now, but I'm still here to help guide your learning journey. What would you like to explore?",
              basic: `I'd love to help you learn about ${responseQuery}! I'm experiencing some connectivity issues with my knowledge base right now, but I can still guide you through this topic. What specific aspect interests you most?`,
              reasoning: `That's a great topic to think through step by step! While I'm having some connectivity issues, I can still help you break down ${responseQuery} into logical components. What's your current understanding of this topic?`,
              research: `Excellent choice for deep research! I'm experiencing some connectivity issues with my research databases, but I can still help you explore ${responseQuery} systematically. What angle would you like to investigate first?`
            };
            
            response = {
              id: 'fallback_' + Date.now(),
              model: 'fallback',
              content: fallbackResponses[determinedMode] || fallbackResponses.basic,
              fullContent: fallbackResponses[determinedMode] || fallbackResponses.basic,
              citations: [],
              metadata: { tokensUsed: 0, responseTime: 0 }
            };
          }

          // Add to conversation history
          updateConversationContext(socket.id, command, response.fullContent || response.content, response.content);
          
          // Create spatial learning panels based on response type
          const panelData = sessionState.isFirstInteraction ? [] : createSpatialPanels(response, determinedMode, responseQuery, model);
          
          // Send response to frontend
          console.log('🔥 DEBUGGING TTS MESSAGE:', {
            messageContent: response.content,
            messageLength: response.content?.length,
            fullContent: response.fullContent,
            fullContentLength: response.fullContent?.length,
            messagePreview: response.content?.substring(0, 100)
          });
          
          socket.emit('voice_response', {
            message: response.content, // Short TTS message
            data: {
              panels: panelData,
              mode: determinedMode,
              model: model,
              query: responseQuery,
              citations: response.citations || [],
              reasoning: response.reasoning || [],
              metadata: response.metadata
            }
          });

          logger.info(`Monday response sent to ${socket.id}:`, {
            service: 'monday-backend',
            mode: determinedMode,
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
function createSpatialPanels(response: any, mode: string, query: string, model: string): any[] {
  console.log('🎯 createSpatialPanels: Creating static panels for response:', {
    mode,
    model,
    query: query.substring(0, 50),
    hasMetadata: !!response.metadata,
    isThinking: response.metadata?.isThinking,
    isResearching: response.metadata?.isResearching,
    hasFullContent: !!response.fullContent,
    fullContentLength: response.fullContent?.length || 0
  })
  
  const panels: any[] = [];
  
  // Handle thinking/researching responses differently
  if (response.metadata?.isThinking || response.metadata?.isResearching) {
    const isThinking = response.metadata?.isThinking;
    const isResearching = response.metadata?.isResearching;
    const actionWord = isThinking ? 'thinking' : 'researching';
    const cleanQuery = query.replace(/^(please\s+)?(think\s+(through\s+|about\s+)?|research\s+|investigate\s+)/i, '').trim();
    
    console.log('🎯 createSpatialPanels: Creating progressive panels for', actionWord, 'process')
    
    // Main panel shows the final research content
    const mainPanel = {
      id: `panel_${Date.now()}_main`,
      type: 'content',
      position: [0, 1.5, -2],
      rotation: [0, 0, 0],
      title: isThinking ? 'Monday\'s Reasoning Process' : 'Monday\'s Research Analysis',
      content: response.content,
      fullContent: response.fullContent,
      isActive: true,
      opacity: 1,
      createdAt: Date.now(),
      model: model,
      isThinking: isThinking,
      isResearching: isResearching,
      citations: response.citations || []
    }
    panels.push(mainPanel)
    console.log('🎯 createSpatialPanels: Created main panel:', mainPanel.id, mainPanel.type)
    
    // Only create extra panels for deep research queries
    if (isResearching) {
      // Thinking process panel on the left
      if (response.thinkingProcess) {
        const thinkingPanel = {
          id: `panel_${Date.now()}_thinking`,
          type: 'content',
          position: [-2.5, 1.5, -1.5],
          rotation: [0, 15, 0],
          title: 'Research Process',
          content: response.thinkingProcess,
          fullContent: response.thinkingProcess,
          reasoning: response.reasoning || [],
          sources: response.sources || [],
          citations: response.citations || [],
          isActive: true,
          opacity: 0.9,
          createdAt: Date.now(),
          model: model
        }
        panels.push(thinkingPanel)
        console.log('🎯 createSpatialPanels: Created thinking process panel:', thinkingPanel.id)
      }
      
      // Progress panel on the right
      const progressPanel = {
        id: `panel_${Date.now()}_progress`,
        type: 'content',
        position: [2.5, 1.5, -1.5],
        rotation: [0, -15, 0],
        title: 'Research Progress',
        content: `I'm researching through ${cleanQuery} step by step. Please wait while I work through this systematically.`,
        fullContent: response.fullContent,
        reasoning: response.reasoning || [],
        sources: response.sources || [],
        citations: response.citations || [],
        isActive: true,
        opacity: 0.9,
        createdAt: Date.now(),
        model: model
      }
      panels.push(progressPanel)
      console.log('🎯 createSpatialPanels: Created progress panel:', progressPanel.id)
    }
    
    console.log('🎯 createSpatialPanels: Total panels created for progressive response:', panels.length)
    return panels;
  }
  
  console.log('🎯 createSpatialPanels: Creating regular panels (not progressive)')
  
  // Regular response handling (non-thinking/researching)
  const panelContent = response.fullContent || response.content;
  
  const mainPanel = {
    id: `panel_${Date.now()}_main`,
    type: 'content',
    position: [-1.5, 1.6, -2],
    rotation: [0, 0, 0],
    title: mode === 'reasoning' ? 'Monday is thinking...' : 
           mode === 'research' ? 'Monday is researching...' :
           mode === 'greeting' ? 'Welcome to Monday' : 
           `Monday: ${query.substring(0, 50)}${query.length > 50 ? '...' : ''}`,
    content: response.content, // Short TTS message
    fullContent: response.fullContent, // Full content for display
    isActive: true,
    opacity: 1,
    createdAt: Date.now(),
    model: model,
    citations: response.citations || []
  };
  panels.push(mainPanel)
  console.log('🎯 createSpatialPanels: Created main panel:', mainPanel.id)
  
  // Reasoning/Research panel (for thinking/researching responses)
  if (response.metadata?.isThinking || response.metadata?.isResearching) {
    const isThinking = response.metadata?.isThinking;
    
    const reasoningPanel = {
      id: `panel_${Date.now()}_reasoning`,
      type: 'reasoning',
      position: [1.5, 1.6, -2],
      rotation: [0, 0, 0],
      title: isThinking ? 'Reasoning Process' : 'Research Analysis',
      content: response.fullContent || 'Analysis in progress...', // Full reasoning/research content
      fullContent: response.fullContent,
      reasoning: response.reasoning || [],
      sources: response.sources || [],
      citations: response.citations || [],
      isActive: false,
      opacity: 1,
      createdAt: Date.now(),
      model: model
    }
    panels.push(reasoningPanel)
    console.log('🎯 createSpatialPanels: Created reasoning panel:', reasoningPanel.id)
  }
  
  console.log('🎯 createSpatialPanels: Total panels created:', panels.length)
  return panels;
}