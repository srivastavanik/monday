const express = require('express');
const { createServer } = require('http');
const { Server } = require('socket.io');
const cors = require('cors');

// Set environment variables directly since .env isn't being read properly
if (!process.env.PERPLEXITY_API_KEY) {
  process.env.PERPLEXITY_API_KEY = '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';
}

// YouTube API key
if (!process.env.YOUTUBE_API_KEY) {
  process.env.YOUTUBE_API_KEY = 'AIzaSyDn_zCV8AGjkQFufH7RDGkSiXD75-2Q39M';
}

const app = express();
const server = createServer(app);

// More debugging for Socket.IO
console.log('Initializing Socket.IO server...');
const io = new Server(server, {
  cors: {
    origin: "*", // Allow all origins for debugging
    methods: ["GET", "POST", "OPTIONS"],
    credentials: true,
    allowedHeaders: ["Content-Type"]
  },
  transports: ['websocket', 'polling'], // Prefer websocket
  pingTimeout: 60000, // 60 seconds
  pingInterval: 25000, // 25 seconds
  allowEIO3: true // Allow older Socket.IO clients
});

// Debug middleware for Socket.IO
io.engine.on('connection_error', (err) => {
  console.log('Socket.IO engine connection error:', err.req);
  console.log('Socket.IO engine error message:', err.message);
  console.log('Socket.IO engine error type:', err.type);
  console.log('Socket.IO engine error context:', err.context);
});

const PORT = 3001;

// Middleware
app.use(cors({
  origin: "*", // Allow all origins for debugging
  credentials: true
}));
app.use(express.json());

// Add root endpoint for testing
app.get('/', (req, res) => {
  res.json({
    message: 'Monday backend server is running',
    socketio: 'enabled',
    port: PORT
  });
});

// Health endpoint
app.get('/health', (req, res) => {
  res.json({
    status: 'healthy',
    timestamp: new Date().toISOString(),
    uptime: process.uptime(),
    perplexityApiConfigured: !!process.env.PERPLEXITY_API_KEY,
    socketConnections: io.engine.clientsCount || 0
  });
});

// Socket.IO test endpoint
app.get('/socket-test', (req, res) => {
  res.json({
    socketIO: 'ready',
    engine: io.engine ? 'initialized' : 'not initialized',
    transports: ['polling', 'websocket'],
    cors: {
      origin: '*',
      credentials: true
    },
    activeConnections: io.engine ? io.engine.clientsCount : 0,
    url: 'ws://localhost:3001'
  });
});

// Enhanced mode detection for voice commands
function detectVoiceMode(command) {
  const commandLower = command.toLowerCase().trim();
  
  // Reasoning mode triggers - comprehensive natural language patterns
  const reasoningTriggers = [
    "think about", "analyze", "reasoning", "figure out", "work through", "solve this",
    "explain why", "help me understand", "break down", "think through", "reason about",
    "logic behind", "walk me through", "make sense of", "explain the reasoning",
    "process this", "work this out", "think it through", "analyze this", "help me think",
    "what's the logic", "how does this work", "can you reason", "step by step",
    "logical analysis", "critical thinking", "problem solving", "analytical approach"
  ];
  
  // Deep research mode triggers - comprehensive investigation patterns
  const researchTriggers = [
    "research into", "investigate", "deep dive", "find information about", "look into",
    "study", "explore", "deep research", "comprehensive analysis", "thorough investigation",
    "research this", "dig deeper", "find out about", "learn everything about",
    "comprehensive study", "detailed research", "in-depth analysis", "full investigation",
    "tell me everything", "complete overview", "detailed explanation", "extensive research",
    "academic research", "scholarly analysis", "comprehensive review", "thorough study"
  ];
  
  // Basic mode triggers - quick search and simple queries
  const basicTriggers = [
    "search the web", "find online", "what is", "define", "quick question",
    "search for", "look up", "find", "basic search", "simple query",
    "quick search", "web search", "google", "find me", "search", "who is",
    "when did", "where is", "how much", "how many", "simple answer",
    "quick answer", "basic question", "fast search", "just tell me"
  ];
  
  // Check reasoning triggers first (highest priority for analytical thinking)
  for (const trigger of reasoningTriggers) {
    if (commandLower.includes(trigger)) {
      return { mode: "reasoning", confidence: 0.95 };
    }
  }
  
  // Check research triggers (high priority for comprehensive investigation)
  for (const trigger of researchTriggers) {
    if (commandLower.includes(trigger)) {
      return { mode: "deep-research", confidence: 0.9 };
    }
  }
  
  // Check basic triggers (standard search queries)
  for (const trigger of basicTriggers) {
    if (commandLower.includes(trigger)) {
      return { mode: "basic", confidence: 0.8 };
    }
  }
  
  // Advanced pattern matching for implicit reasoning requests
  if (commandLower.match(/\b(why|how|because|reason|cause|effect|impact|consequence)\b/)) {
    return { mode: "reasoning", confidence: 0.7 };
  }
  
  // Advanced pattern matching for research requests
  if (commandLower.match(/\b(history|background|origin|development|evolution|comprehensive|detailed)\b/)) {
    return { mode: "deep-research", confidence: 0.6 };
  }
  
  // For conversational queries without specific triggers, use reasoning mode for better responses
  if (commandLower.match(/\b(can you|could you|would you|tell me about|explain|describe|discuss)\b/)) {
    return { mode: "reasoning", confidence: 0.5 };
  }
  
  // Default to basic mode for simple queries
  return { mode: "basic", confidence: 0.3 };
}

// Perplexity API integration
async function queryPerplexity(prompt, mode = 'basic') {
  const modelMap = {
    'basic': 'sonar',
    'reasoning': 'sonar-reasoning-pro', 
    'deep-research': 'sonar-deep-research'
  };

  const model = modelMap[mode] || modelMap['basic'];
  
  console.log(`Querying Perplexity API with model: ${model} for prompt: "${prompt.substring(0, 50)}..."`);
  
  try {
    const response = await fetch('https://api.perplexity.ai/chat/completions', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${process.env.PERPLEXITY_API_KEY}`,
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        model: model,
        messages: [
          {
            role: 'system',
            content: `You are Monday, an advanced AI learning companion and educational assistant designed for voice interactions. You have a friendly, knowledgeable personality and specialize in making complex topics easy to understand through natural conversation.

PERSONALITY:
- Enthusiastic about learning and teaching
- Patient and encouraging  
- Conversational and natural (like talking to a knowledgeable friend)
- Uses clear, educational explanations
- Provides practical examples and analogies
- Encourages curiosity and deeper exploration

RESPONSE STYLE FOR VOICE INTERACTION:
- Respond naturally as if speaking directly to the user
- Keep responses conversational but informative
- Use "you" and "I" appropriately for natural dialogue
- Break down complex concepts into digestible parts
- Use analogies and real-world examples
- Suggest follow-up questions or related topics to explore
- Keep responses engaging but not overly long for voice

MODE-SPECIFIC BEHAVIOR:
${mode === 'reasoning' ? `- REASONING MODE: Provide thoughtful analysis and logical explanations
- Walk through your thinking process step by step
- Explain the "why" behind concepts and phenomena
- Use analytical frameworks and reasoning approaches
- Connect ideas and show relationships between concepts` : ''}
${mode === 'deep-research' ? `- DEEP RESEARCH MODE: Provide comprehensive, well-researched information
- Draw from multiple perspectives and sources
- Provide thorough background and context
- Include relevant details and nuances
- Offer deeper insights and connections to broader topics` : ''}
${mode === 'basic' ? `- BASIC MODE: Provide clear, concise answers and explanations
- Focus on the most important and relevant information
- Keep explanations accessible and easy to understand
- Be direct while still being conversational and helpful` : ''}

CONVERSATION CONTEXT:
- This is part of an ongoing voice conversation
- The user just asked you a question directly
- Respond as Monday, the helpful learning companion
- Make your response feel like a natural part of the conversation

Always respond as Monday in a conversational, voice-appropriate manner. You are having a real-time conversation with the user.`
          },
          {
            role: 'user',
            content: prompt
          }
        ],
        max_tokens: mode === 'deep-research' ? 2048 : 1024,
        temperature: 0.2,
        top_p: 0.9,
        stream: false
      })
    });

    if (!response.ok) {
      const errorText = await response.text();
      console.error(`Perplexity API error: ${response.status} - ${errorText}`);
      throw new Error(`Perplexity API error: ${response.status} - ${errorText}`);
    }

    const data = await response.json();
    console.log('Perplexity API response received successfully');
    
    // Filter out thinking tags from the response
    let content = data.choices[0]?.message?.content || 'I apologize, but I didn\'t receive a proper response. Could you please try asking your question again?';
    content = filterThinkingTags(content);
    
    return {
      content: content,
      model: model,
      usage: data.usage
    };
  } catch (error) {
    console.error('Perplexity API error:', error);
    throw error; // Re-throw to handle in the calling function
  }
}

// Function to filter out thinking tags and internal reasoning
function filterThinkingTags(content) {
  if (!content) return content;
  
  // Remove <think>...</think> blocks completely
  content = content.replace(/<think>[\s\S]*?<\/think>/gi, '');
  
  // Remove any remaining thinking artifacts
  content = content.replace(/^\s*thinking:[\s\S]*?(?=\n\n|\n[A-Z]|$)/im, '');
  content = content.replace(/^\s*analysis:[\s\S]*?(?=\n\n|\n[A-Z]|$)/im, '');
  
  // Clean up extra whitespace
  content = content.replace(/\n\s*\n\s*\n/g, '\n\n');
  content = content.trim();
  
  return content;
}

// Real YouTube Data API search function
async function searchEducationalVideos(query, mode) {
  console.log(`Searching YouTube API for educational videos with query: "${query}" in ${mode} mode`);
  
  if (!process.env.YOUTUBE_API_KEY) {
    console.error('YouTube API key not configured');
    return null;
  }
  
  try {
    // Create educational search terms based on the query and mode
    const educationalTerms = ['tutorial', 'explained', 'lesson', 'course', 'learn', 'education'];
    const searchQuery = `${query} ${educationalTerms[Math.floor(Math.random() * educationalTerms.length)]}`;
    
    // Call YouTube Data API v3
    const apiUrl = `https://www.googleapis.com/youtube/v3/search`;
    const params = new URLSearchParams({
      part: 'snippet',
      q: searchQuery,
      key: process.env.YOUTUBE_API_KEY,
      type: 'video',
      order: 'relevance',
      maxResults: '3', // Get 3 videos (1 main + 2 related)
      safeSearch: 'moderate',
      videoEmbeddable: 'true',
      videoDuration: 'medium', // 4-20 minutes are good for educational content
      relevanceLanguage: 'en',
      regionCode: 'US'
    });

    const response = await fetch(`${apiUrl}?${params.toString()}`);
    
    if (!response.ok) {
      console.error(`YouTube API error: ${response.status} - ${response.statusText}`);
      throw new Error(`YouTube API request failed: ${response.status}`);
    }

    const data = await response.json();
    
    if (!data.items || data.items.length === 0) {
      console.log('No videos found, falling back to default');
      return getFallbackVideo();
    }

    // Get the first video as the main result
    const mainVideo = data.items[0];
    
    console.log(`Found YouTube video: "${mainVideo.snippet.title}" by ${mainVideo.snippet.channelTitle}`);
    
    return {
      id: mainVideo.id.videoId,
      title: mainVideo.snippet.title,
      channel: mainVideo.snippet.channelTitle,
      description: mainVideo.snippet.description,
      thumbnail: mainVideo.snippet.thumbnails.medium?.url || mainVideo.snippet.thumbnails.default?.url,
      publishedAt: new Date(mainVideo.snippet.publishedAt).toLocaleDateString(),
      // Also return related videos for the frontend
      relatedVideos: data.items.slice(1).map(video => ({
        id: video.id.videoId,
        title: video.snippet.title,
        channel: video.snippet.channelTitle,
        description: video.snippet.description,
        thumbnail: video.snippet.thumbnails.medium?.url || video.snippet.thumbnails.default?.url,
        publishedAt: new Date(video.snippet.publishedAt).toLocaleDateString()
      }))
    };

  } catch (error) {
    console.error('YouTube API search failed:', error);
    
    // Try a simpler search without educational terms
    try {
      console.log('Retrying with simpler search...');
      const simpleApiUrl = `https://www.googleapis.com/youtube/v3/search`;
      const simpleParams = new URLSearchParams({
        part: 'snippet',
        q: query,
        key: process.env.YOUTUBE_API_KEY,
        type: 'video',
        order: 'relevance',
        maxResults: '1',
        safeSearch: 'moderate',
        videoEmbeddable: 'true'
      });

      const simpleResponse = await fetch(`${simpleApiUrl}?${simpleParams.toString()}`);
      
      if (simpleResponse.ok) {
        const simpleData = await simpleResponse.json();
        if (simpleData.items && simpleData.items.length > 0) {
          const video = simpleData.items[0];
          console.log(`Fallback search found: "${video.snippet.title}"`);
          return {
            id: video.id.videoId,
            title: video.snippet.title,
            channel: video.snippet.channelTitle,
            description: video.snippet.description,
            thumbnail: video.snippet.thumbnails.medium?.url || video.snippet.thumbnails.default?.url,
            publishedAt: new Date(video.snippet.publishedAt).toLocaleDateString(),
            relatedVideos: []
          };
        }
      }
    } catch (retryError) {
      console.error('Retry search also failed:', retryError);
    }
    
    // Final fallback to a known good educational video
    return getFallbackVideo();
  }
}

// Fallback video when YouTube API fails
function getFallbackVideo() {
  console.log('Using fallback educational video');
  return {
    id: "aircAruvnKk", // 3Blue1Brown Neural Networks - very popular and reliable
    title: "But what is a Neural Network? | Deep learning, chapter 1",
    channel: "3Blue1Brown",
    description: "Subscribe for new videos every Friday! Neural networks can seem like a black box...",
    thumbnail: "https://img.youtube.com/vi/aircAruvnKk/mqdefault.jpg",
    publishedAt: "Oct 5, 2017",
    relatedVideos: [
      {
        id: "IHZwWFHWa-w",
        title: "Gradient descent, how neural networks learn | Deep learning, chapter 2",
        channel: "3Blue1Brown",
        description: "How do neural networks learn? In this video, we dive into the gradient descent algorithm...",
        thumbnail: "https://img.youtube.com/vi/IHZwWFHWa-w/mqdefault.jpg",
        publishedAt: "Oct 16, 2017"
      },
      {
        id: "Ilg3gGewQ5U", 
        title: "What is backpropagation really doing? | Deep learning, chapter 3",
        channel: "3Blue1Brown",
        description: "Backpropagation is the heart of neural network training...",
        thumbnail: "https://img.youtube.com/vi/Ilg3gGewQ5U/mqdefault.jpg",
        publishedAt: "Nov 3, 2017"
      }
    ]
  };
}

// Socket.IO connection handling
io.on('connection', (socket) => {
  console.log('=== NEW CLIENT CONNECTION ===');
  console.log('Client connected:', socket.id);
  console.log('Client address:', socket.handshake.address);
  console.log('Transport:', socket.conn.transport.name);
  console.log('Total clients:', io.engine.clientsCount);
  console.log('============================');

  socket.on('voice_command', async (data) => {
    console.log('Voice command received:', data);
    
    const { command, conversationActive, isExplicitTrigger, isActivation } = data;
    
    try {
      // If this is just the activation (Hey Monday without additional content), respond briefly
      if (isActivation && (!command || command.trim().length === 0)) {
        socket.emit('voice_response', {
          message: "Hi! I'm Monday, your learning assistant powered by Perplexity Sonar. What would you like to dive into today?",
          data: {
            panels: [],
            mode: 'basic',
            model: 'sonar',
            query: 'activation',
            citations: [],
            reasoning: [],
            metadata: {
              tokensUsed: 0,
              responseTime: Date.now(),
              isActivation: true
            }
          }
        });
        return;
      }

      // For all other cases (including Hey Monday + command), process with Perplexity
      const modeDetection = detectVoiceMode(command);
      const mode = modeDetection.mode;
      
      console.log(`Processing command in ${mode} mode (confidence: ${modeDetection.confidence}): "${command}"`);
      
      // Query Perplexity API with the actual command (removed progress update delay)
      const perplexityResponse = await queryPerplexity(command, mode);
      
      // Get relevant YouTube video
      const youtubeVideo = await searchEducationalVideos(command, mode);
      
      console.log(`Generated response: "${perplexityResponse.content.substring(0, 100)}..."`);
      
      // Create spatial panels
      const panels = [{
        id: `panel_${Date.now()}_main`,
        type: 'content',
        position: [0, 1.5, -2],
        rotation: [0, 0, 0],
        title: `Monday: ${command.substring(0, 50)}`,
        content: perplexityResponse.content,
        fullContent: perplexityResponse.content,
        isActive: true,
        opacity: 1,
        createdAt: Date.now(),
        model: perplexityResponse.model,
        citations: []
      }];
      
      // Send response
      socket.emit('voice_response', {
        message: perplexityResponse.content,
        data: {
          panels: panels,
          mode: mode,
          model: perplexityResponse.model,
          query: command,
          citations: [],
          reasoning: [],
          metadata: {
            youtubeVideoId: youtubeVideo?.id,
            youtubeTitle: youtubeVideo?.title,
            youtubeChannel: youtubeVideo?.channel,
            youtubeDescription: youtubeVideo?.description,
            youtubeThumbnail: youtubeVideo?.thumbnail,
            youtubePublishedAt: youtubeVideo?.publishedAt,
            youtubeRelatedVideos: youtubeVideo?.relatedVideos || [],
            tokensUsed: perplexityResponse.usage?.total_tokens || 0,
            responseTime: Date.now()
          }
        }
      });
      
    } catch (error) {
      console.error('Error processing voice command:', error);
      socket.emit('voice_error', {
        error: `I encountered an issue: ${error.message}. Please try again.`
      });
    }
  });

  // Heartbeat handler to keep connection alive
  socket.on('heartbeat', (data) => {
    socket.emit('heartbeat_ack', { timestamp: Date.now(), received: data.timestamp });
  });

  socket.on('disconnect', () => {
    console.log('Client disconnected:', socket.id);
  });
});

server.listen(PORT, () => {
  console.log(`Monday backend server running on port ${PORT}`);
  console.log(`Health endpoint: http://localhost:${PORT}/health`);
  console.log('WebSocket server ready for connections');
  console.log(`Perplexity API configured: ${!!process.env.PERPLEXITY_API_KEY}`);
  console.log(`Perplexity API key (first 10 chars): ${process.env.PERPLEXITY_API_KEY?.substring(0, 10)}...`);
}); 