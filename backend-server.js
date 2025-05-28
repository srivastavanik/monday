const express = require('express');
const { createServer } = require('http');
const { Server } = require('socket.io');
const cors = require('cors');

// Set environment variables directly since .env isn't being read properly
if (!process.env.PERPLEXITY_API_KEY) {
  process.env.PERPLEXITY_API_KEY = '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';
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

// Enhanced Perplexity API integration with streaming support
async function queryPerplexityStreaming(prompt, mode = 'basic', socket = null) {
  const modelMap = {
    'basic': 'sonar',
    'reasoning': 'sonar-reasoning-pro', 
    'deep-research': 'sonar-deep-research'
  };

  const model = modelMap[mode] || modelMap['basic'];
  
  console.log(`Querying Perplexity API with model: ${model} for prompt: "${prompt.substring(0, 50)}..."`);
  
  // Emit initial progress for reasoning/research modes
  if (socket && (mode === 'reasoning' || mode === 'deep-research')) {
    socket.emit(mode === 'reasoning' ? 'reasoning_progress' : 'research_progress', {
      type: mode === 'reasoning' ? 'reasoning_progress' : 'research_progress',
      update: mode === 'reasoning' ? 'Initializing reasoning engine...' : 'Starting comprehensive research...',
      data: {
        sources: [],
        reasoning: [],
        citations: []
      }
    });
  }
  
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
        max_tokens: mode === 'deep-research' ? 4096 : mode === 'reasoning' ? 2048 : 1024,
        temperature: 0.2,
        top_p: 0.9,
        stream: mode !== 'basic', // Stream for reasoning and research modes
        return_citations: mode === 'deep-research', // Request citations for deep research
        return_sources: mode === 'deep-research', // Request sources for deep research
        search_domain_filter: mode === 'deep-research' ? [] : undefined, // No domain restrictions for deep research
        search_recency_filter: mode === 'deep-research' ? 'month' : undefined // Recent sources for deep research
      })
    });

    if (!response.ok) {
      const errorText = await response.text();
      console.error(`Perplexity API error: ${response.status} - ${errorText}`);
      throw new Error(`Perplexity API error: ${response.status} - ${errorText}`);
    }

    // Handle streaming response for reasoning/research modes
    if (mode !== 'basic' && response.body) {
      console.log('Processing streaming response...');
      
      // Read the full response
      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let fullContent = '';
      let sources = [];
      let citations = [];
      let reasoning = [];
      let lastProgressUpdate = '';
      let seenSources = new Set();
      let seenCitations = new Set();
      let thinkingSteps = [];
      let currentThought = '';
      
      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        
        const chunk = decoder.decode(value);
        const lines = chunk.split('\n');
        
        for (const line of lines) {
          if (line.startsWith('data: ')) {
            const data = line.slice(6);
            if (data === '[DONE]') continue;
            
            try {
              const parsed = JSON.parse(data);
              
              // Extract content
              if (parsed.choices?.[0]?.delta?.content) {
                const deltaContent = parsed.choices[0].delta.content;
                fullContent += deltaContent;
                currentThought += deltaContent;
                
                // For reasoning mode, extract thinking steps from the content
                if (mode === 'reasoning' && deltaContent.includes('\n')) {
                  const thoughtLines = currentThought.split('\n').filter(l => l.trim());
                  for (const thought of thoughtLines) {
                    if (thought.trim() && !thinkingSteps.includes(thought.trim())) {
                      thinkingSteps.push(thought.trim());
                      
                      // Emit reasoning step
                      if (socket) {
                        socket.emit('reasoning_progress', {
                          type: 'reasoning_progress',
                          update: thought.trim(),
                          data: {
                            reasoning: thinkingSteps.map((t, i) => ({
                              content: t,
                              step: i + 1
                            }))
                          }
                        });
                      }
                    }
                  }
                  currentThought = '';
                }
              }
              
              // Extract sources from metadata (for deep research)
              if (parsed.choices?.[0]?.delta?.search_results) {
                const searchResults = parsed.choices[0].delta.search_results;
                for (const result of searchResults) {
                  const sourceKey = result.url || result.title;
                  if (sourceKey && !seenSources.has(sourceKey)) {
                    seenSources.add(sourceKey);
                    const source = {
                      title: result.title || 'Research Source',
                      url: result.url || '',
                      description: result.snippet || result.description || '',
                      domain: result.domain || new URL(result.url || 'http://example.com').hostname
                    };
                    sources.push(source);
                    
                    // Emit source update
                    if (socket && mode === 'deep-research') {
                      socket.emit('research_progress', {
                        type: 'research_progress',
                        update: `Analyzing: ${source.domain}`,
                        data: {
                          sources: sources,
                          citations: citations
                        }
                      });
                    }
                  }
                }
              }
              
              // Extract web results (alternative source format)
              if (parsed.choices?.[0]?.delta?.web_results) {
                const webResults = parsed.choices[0].delta.web_results;
                for (const result of webResults) {
                  const sourceKey = result.url || result.title;
                  if (sourceKey && !seenSources.has(sourceKey)) {
                    seenSources.add(sourceKey);
                    const source = {
                      title: result.title || 'Web Source',
                      url: result.url || '',
                      description: result.snippet || '',
                      domain: new URL(result.url || 'http://example.com').hostname
                    };
                    sources.push(source);
                    
                    // Emit source update
                    if (socket && mode === 'deep-research') {
                      socket.emit('research_progress', {
                        type: 'research_progress',
                        update: `Researching: ${source.domain}`,
                        data: {
                          sources: sources,
                          citations: citations
                        }
                      });
                    }
                  }
                }
              }
              
              // Extract citations
              if (parsed.choices?.[0]?.delta?.citations) {
                for (const citation of parsed.choices[0].delta.citations) {
                  const citationKey = citation.url || citation.title;
                  if (citationKey && !seenCitations.has(citationKey)) {
                    seenCitations.add(citationKey);
                    citations.push({
                      title: citation.title || citation.text || 'Citation',
                      url: citation.url || '',
                      author: citation.author || '',
                      year: citation.year || ''
                    });
                  }
                }
              }
              
            } catch (e) {
              // Ignore parse errors for incomplete chunks
              console.log('Parse error (expected for partial chunks):', e.message);
            }
          }
        }
      }
      
      // Process any remaining thought
      if (currentThought.trim() && mode === 'reasoning' && !thinkingSteps.includes(currentThought.trim())) {
        thinkingSteps.push(currentThought.trim());
        if (socket) {
          socket.emit('reasoning_progress', {
            type: 'reasoning_progress',
            update: currentThought.trim(),
            data: {
              reasoning: thinkingSteps.map((t, i) => ({
                content: t,
                step: i + 1
              }))
            }
          });
        }
      }
      
      // Filter out thinking tags from the response
      fullContent = filterThinkingTags(fullContent);
      
      // Final progress update
      if (socket) {
        const finalUpdate = mode === 'reasoning' 
          ? 'Reasoning complete. Formulating response...'
          : `Research complete. Analyzed ${sources.length} sources.`;
          
        socket.emit(mode === 'reasoning' ? 'reasoning_progress' : 'research_progress', {
          type: mode === 'reasoning' ? 'reasoning_progress' : 'research_progress',
          update: finalUpdate,
          data: {
            sources: sources,
            citations: citations,
            reasoning: mode === 'reasoning' ? thinkingSteps.map((t, i) => ({
              content: t,
              step: i + 1
            })) : []
          }
        });
      }
      
      return {
        content: fullContent || 'I apologize, but I didn\'t receive a proper response. Could you please try asking your question again?',
        model: model,
        usage: { total_tokens: 0 }, // Streaming doesn't provide token usage
        sources: sources,
        citations: citations
      };
    } else {
      // Non-streaming response for basic mode
      const data = await response.json();
      console.log('Perplexity API response received successfully');
      
      // Filter out thinking tags from the response
      let content = data.choices[0]?.message?.content || 'I apologize, but I didn\'t receive a proper response. Could you please try asking your question again?';
      content = filterThinkingTags(content);
      
      return {
        content: content,
        model: model,
        usage: data.usage,
        sources: [],
        citations: []
      };
    }
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

// YouTube search function
function searchEducationalVideos(query, mode) {
  // Educational video suggestions based on query
  const videoSuggestions = [
    { id: "8hly31xKli0", title: "Binary Search Algorithm - Explained Simply", channel: "CS Dojo" },
    { id: "D6xkbGLQesk", title: "Binary Search Tree Implementation", channel: "GeeksforGeeks" },
    { id: "5xlIPT1FRcA", title: "Data Structures: Binary Search Trees", channel: "MIT OpenCourseWare" },
    { id: "pYT9F8_LFTM", title: "Algorithms: Binary Search", channel: "Khan Academy Computer Science" },
    { id: "P3YID7liBug", title: "Binary Search Tree Visualization", channel: "Coursera" },
    { id: "kPRA0W1kECg", title: "Introduction to Machine Learning", channel: "MIT OpenCourseWare" },
    { id: "aircAruvnKk", title: "Neural Networks Explained", channel: "3Blue1Brown" },
    { id: "WUvTyaaNkzM", title: "Calculus Fundamentals", channel: "Khan Academy" },
    { id: "fNk_zzaMoSs", title: "Linear Algebra Essence", channel: "3Blue1Brown" },
    { id: "YAXLy4jNhzI", title: "Physics Concepts Visualized", channel: "MinutePhysics" }
  ];

  // More intelligent video matching
  const queryLower = query.toLowerCase();
  const queryWords = queryLower.split(/\s+/);
  
  // Find videos that match query keywords
  let relevantVideos = videoSuggestions.filter(video => {
    const videoTitle = video.title.toLowerCase();
    return queryWords.some(word => 
      word.length > 2 && videoTitle.includes(word)
    );
  });

  // If no direct matches, find videos with related topics
  if (relevantVideos.length === 0) {
    const topicKeywords = {
      'binary': ['binary', 'search', 'tree', 'algorithm'],
      'machine': ['machine', 'learning', 'neural', 'ai'],
      'math': ['calculus', 'algebra', 'math'],
      'physics': ['physics', 'science'],
      'programming': ['algorithm', 'code', 'programming']
    };

    for (const [topic, keywords] of Object.entries(topicKeywords)) {
      if (queryWords.some(word => keywords.includes(word))) {
        relevantVideos = videoSuggestions.filter(video =>
          keywords.some(keyword => video.title.toLowerCase().includes(keyword))
        );
        break;
      }
    }
  }

  return relevantVideos.length > 0 ? relevantVideos[0] : videoSuggestions[0];
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
      
      // Query Perplexity API with streaming support for reasoning/research modes
      const perplexityResponse = await queryPerplexityStreaming(command, mode, socket);
      
      // Get relevant YouTube video
      const youtubeVideo = searchEducationalVideos(command, mode);
      
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
          citations: perplexityResponse.citations || [],
          reasoning: [],
          sources: perplexityResponse.sources || [],
          metadata: {
            youtubeVideoId: youtubeVideo.id,
            youtubeTitle: youtubeVideo.title,
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