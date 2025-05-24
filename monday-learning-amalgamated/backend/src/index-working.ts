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
    origin: process.env.FRONTEND_URL || "http://localhost:3000",
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

// Basic middleware
app.use(cors({
  origin: process.env.FRONTEND_URL || "http://localhost:3000",
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

// Helper functions
function detectQueryMode(command: string): 'greeting' | 'basic' | 'reasoning' | 'research' {
  const lowerCommand = command.toLowerCase();
  
  if (lowerCommand.includes('hello') || lowerCommand.includes('hi') || lowerCommand.includes('hey')) {
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

function getConversationContext(socketId: string): string[] {
  return conversationContexts.get(socketId) || [];
}

function updateConversationContext(socketId: string, userInput: string, response: string) {
  const context = getConversationContext(socketId);
  context.push(`User: ${userInput}`);
  context.push(`Monday: ${response}`);
  
  if (context.length > 20) {
    context.splice(0, context.length - 20);
  }
  
  conversationContexts.set(socketId, context);
}

function createSpatialPanels(response: any, mode: string, query: string): any[] {
  const panels: any[] = [];
  
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

// Socket.IO connection handling
io.on('connection', (socket) => {
  logger.info(`Client connected: ${socket.id}`);
  
  conversationContexts.set(socket.id, []);
  
  socket.on('voice_command', async (data) => {
    try {
      const { command } = data;
      logger.info(`Voice command received: "${command}" from ${socket.id}`);
      
      const normalizedCommand = command.toLowerCase().trim();
      const sessionState = sessionStates.get(socket.id) || { 
        isInConversation: false, 
        lastActivity: Date.now() 
      };
      
      const isMondayActivation = normalizedCommand.includes('monday');
      const isInActiveConversation = sessionState.isInConversation && 
        (Date.now() - sessionState.lastActivity < 300000);
      
      if (isMondayActivation || isInActiveConversation) {
        const query = command.replace(/monday,?\s*/i, '').trim();
        const mode = detectQueryMode(command);
        const context = getConversationContext(socket.id);
        
        sessionStates.set(socket.id, {
          isInConversation: true,
          lastActivity: Date.now(),
          currentTopic: query || sessionState.currentTopic
        });
        
        try {
          let response;
          
          if (perplexityService) {
            // Use Perplexity if available
            if (mode === 'greeting' && query.length < 10) {
              response = await perplexityService.processQuery({
                query: "The user just said hello. Greet them warmly as Monday, introduce yourself as their AI learning companion, and ask what they'd like to explore or learn about today. Mention briefly that you are powered by Perplexity. Keep it conversational and under 3 sentences.",
                mode: 'basic',
                context: context.slice(-4),
                sessionId: socket.id
              });
            } else if (query.length > 0) {
              const fullQuery = context.length > 0 
                ? `Previous conversation context: ${context.slice(-4).join(' ')}\n\nCurrent user question: ${query}`
                : query;
                
              response = await perplexityService.processQuery({
                query: fullQuery,
                mode: mode === 'basic' ? 'basic' : mode === 'reasoning' ? 'reasoning' : 'research',
                context: context.slice(-6),
                sessionId: socket.id
              });
            } else {
              response = await perplexityService.processQuery({
                query: "The user said 'Monday' but didn't ask anything specific. Politely ask what they'd like to learn about or explore. Be encouraging and suggest some interesting topics.",
                mode: 'basic',
                context: context.slice(-2),
                sessionId: socket.id
              });
            }
          } else {
            // Fallback responses if Perplexity isn't available
            const fallbackResponses = {
              greeting: "Hello! I'm Monday, your AI learning companion. I'd love to help you explore any topic you're curious about. What would you like to learn today?",
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
          
          updateConversationContext(socket.id, command, response.content);
          const panelData = createSpatialPanels(response, mode, query);
          
          const responseData = {
            type: mode === 'greeting' ? 'greeting' : 'learning_response',
            message: response.content,
            action: mode === 'research' ? 'show_research_panel' : mode === 'reasoning' ? 'show_reasoning_panel' : 'show_info_panel',
            data: {
              title: mode === 'greeting' ? 'Welcome to Monday' : `Learning: ${query}`,
              content: response.content,
              citations: response.citations || [],
              reasoning: response.reasoning || [],
              sources: response.sources || [],
              mode: mode,
              timestamp: new Date().toISOString(),
              metadata: response.metadata,
              panels: panelData,
              conversationActive: true
            }
          };
          
          socket.emit('voice_response', responseData);
          logger.info(`Monday response sent to ${socket.id}`, {
            mode,
            responseLength: response.content.length,
            citationCount: response.citations?.length || 0,
            tokensUsed: response.metadata.tokensUsed,
            panelCount: panelData.length
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
        logger.info(`Ignoring command without Monday activation: "${command}"`);
      }
      
    } catch (error) {
      logger.error(`Voice command error for ${socket.id}:`, error);
      socket.emit('voice_error', { 
        error: 'Failed to process voice command',
        timestamp: new Date().toISOString()
      });
    }
  });
  
  socket.on('spatial_command', (data) => {
    logger.info(`Spatial command from ${socket.id}:`, data);
    socket.emit('spatial_response', { 
      success: true, 
      message: 'Spatial command processed',
      timestamp: new Date().toISOString()
    });
  });
  
  socket.on('session_start', (data) => {
    logger.info(`Session started for ${socket.id}:`, data);
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

process.on('uncaughtException', (error) => {
  logger.error('Uncaught Exception:', error);
  process.exit(1);
});

process.on('unhandledRejection', (reason, promise) => {
  logger.error('Unhandled Rejection at:', promise, 'reason:', reason);
  process.exit(1);
});

console.log('[BACKEND LOG] All setup completed - ready for connections!'); 