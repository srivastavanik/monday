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

// Test Perplexity API endpoint
app.get('/test-perplexity', async (req, res) => {
  console.log('Testing Perplexity API...');
  try {
    const testPrompt = 'Say hello and tell me you are Monday, the AI learning assistant.';
    const response = await queryPerplexityStreaming(testPrompt, 'basic', null);
    console.log('Test response:', response);
    res.json({
      success: true,
      response: response.content,
      model: response.model
    });
  } catch (error) {
    console.error('Perplexity API test failed:', error);
    res.status(500).json({
      success: false,
      error: error.message,
      stack: error.stack
    });
  }
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

  // Enhanced token limits for better performance
  const tokenLimits = {
    'basic': 2000,
    'reasoning': 8000, // Increased for better reasoning
    'deep-research': 12000 // Increased for comprehensive research
  };

  const model = modelMap[mode] || modelMap['basic'];
  const maxTokens = tokenLimits[mode] || tokenLimits['basic'];
  
  console.log(`Querying Perplexity API with model: ${model}, tokens: ${maxTokens} for prompt: "${prompt.substring(0, 50)}..."`);
  
  const shouldStream = mode !== 'basic';
  
  // For reasoning and deep-research modes, immediately send an initial progress update
  if (socket && shouldStream) {
    socket.emit(mode === 'reasoning' ? 'reasoning_progress' : 'research_progress', {
      type: mode === 'reasoning' ? 'reasoning_progress' : 'research_progress',
      update: mode === 'reasoning' ? 'Starting analytical reasoning process...' : 'Initiating comprehensive research...',
      data: { reasoning: [], sources: [], citations: [] }
    });
  }

  const timeoutId = setTimeout(() => {
    if (abortController) {
      abortController.abort();
    }
  }, 120000); // 2 minute timeout

  const abortController = new AbortController();

  try {
    const response = await fetch('https://api.perplexity.ai/chat/completions', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${process.env.PERPLEXITY_API_KEY}`,
        'Content-Type': 'application/json',
        'Accept': shouldStream ? 'text/event-stream' : 'application/json'
      },
      signal: abortController.signal,
      body: JSON.stringify({
        model: model,
        messages: [
          {
            role: 'system',
            content: `You are Monday, an intelligent learning assistant powered by Perplexity Sonar. You provide helpful, accurate, and engaging responses in a conversational tone.

MODE-SPECIFIC BEHAVIOR:
${mode === 'reasoning' ? `- REASONING MODE: Provide thoughtful analysis and logical explanations
- Walk through your thinking process step by step
- Begin with a <think> section containing your detailed reasoning
- Explain the "why" behind concepts and phenomena
- Use analytical frameworks and reasoning approaches
- Connect ideas and show relationships between concepts` : ''}
${mode === 'deep-research' ? `- DEEP RESEARCH MODE: Provide comprehensive, well-researched information
- Draw from multiple perspectives and sources
- Document your research process and key findings
- Begin with a <research> section outlining your research approach
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

OUTPUT FORMAT FOR REASONING & RESEARCH:
- For reasoning & deep research modes, ALWAYS include your thinking process
- In reasoning mode, wrap detailed reasoning in <think>...</think> tags
- In deep research mode, wrap research steps in <research>...</research> tags
- This thinking/research section will be displayed to the user separately
- After your thinking/research section, provide your final response

Always respond as Monday in a conversational, voice-appropriate manner. You are having a real-time conversation with the user.`
          },
          {
            role: 'user',
            content: prompt
          }
        ],
        max_tokens: maxTokens,
        temperature: 0.2,
        top_p: 0.9,
        stream: shouldStream,
        return_citations: true,
        return_sources: true
      })
    });

    clearTimeout(timeoutId);
    console.log('Perplexity API response status:', response.status);

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Perplexity API error: ${response.status} - ${errorText}`);
    }

    if (!shouldStream) {
      const data = await response.json();
      let content = data.choices[0]?.message?.content || 'No response content.';
      
      // Extract thinking process for basic mode (if any)
      const thinkingProcess = extractThinkingProcess(content);
      content = filterThinkingTags(content);
      
      return {
        content: content,
        model: model,
        usage: data.usage,
        sources: data.sources || [],
        citations: data.citations || [],
        thinking: thinkingProcess
      };
    }
    
    // Enhanced SSE Streaming handling
    if (shouldStream && response.body) {
      console.log('Processing SSE stream for advanced mode...');
      const reader = response.body.getReader();
      const decoder = new TextDecoder();
      let sseBuffer = "";
      let fullContent = "";
      let thinkingProcess = "";
      let sources = [];
      let citations = [];
      let reasoning = [];
      let progressCount = 0;
      let isCollectingThinking = false;
      let currentThinkingChunk = "";

      while (true) {
        const { done, value } = await reader.read();
        if (done) {
          console.log("Stream finished.");
          break;
        }
        
        sseBuffer += decoder.decode(value, { stream: true });

        let eventEndIndex;
        while ((eventEndIndex = sseBuffer.indexOf('\n\n')) !== -1) {
          const eventString = sseBuffer.substring(0, eventEndIndex);
          sseBuffer = sseBuffer.substring(eventEndIndex + 2);
          const lines = eventString.split('\n');
          for (const line of lines) {
            if (line.startsWith('data: ')) {
              const jsonData = line.substring(6);
              if (jsonData.trim() === '[DONE]') {
                console.log("SSE stream signalled [DONE]");
                continue;
              }
              if (jsonData.trim() === "") continue;

              try {
                const parsedEvent = JSON.parse(jsonData);

                if (parsedEvent.choices && parsedEvent.choices[0] && parsedEvent.choices[0].delta && parsedEvent.choices[0].delta.content) {
                  const deltaContent = parsedEvent.choices[0].delta.content;
                  fullContent += deltaContent;
                  progressCount++;

                  // Enhanced thinking process tracking
                  if (deltaContent.includes("<think>") || deltaContent.includes("<research>")) {
                    isCollectingThinking = true;
                  }
                  
                  if (isCollectingThinking) {
                    currentThinkingChunk += deltaContent;
                    
                    if (deltaContent.includes("</think>") || deltaContent.includes("</research>")) {
                      isCollectingThinking = false;
                      thinkingProcess += currentThinkingChunk;
                      
                      // Extract thinking content
                      const thinkMatch = /<think>([\s\S]*?)<\/think>|<research>([\s\S]*?)<\/research>/.exec(currentThinkingChunk);
                      if (thinkMatch && thinkMatch[1]) {
                        const thinkStep = thinkMatch[1].trim();
                        reasoning.push({ content: thinkStep });
                        
                        // Send real-time thinking/reasoning to client
                        if (socket) {
                          const progressType = mode === 'reasoning' ? 'reasoning_progress' : 'research_progress';
                          socket.emit(progressType, {
                            type: progressType,
                            update: mode === 'reasoning' ? 
                              `Reasoning step ${reasoning.length}: ${thinkStep.substring(0, 100)}...` :
                              `Research insight ${reasoning.length}: ${thinkStep.substring(0, 100)}...`,
                            data: { reasoning, sources, citations }
                          });
                        }
                      }
                      
                      currentThinkingChunk = "";
                    }
                  }

                  // Send periodic progress updates
                  if (socket && progressCount % 3 === 0) {
                     let updateMsg = "Processing response...";
                     if (mode === 'reasoning') updateMsg = `Analyzing: ${deltaContent.substring(0, 50).replace(/<[^>]*>/g, '')}...`;
                     if (mode === 'deep-research') updateMsg = `Researching: ${deltaContent.substring(0, 50).replace(/<[^>]*>/g, '')}...`;
                    
                    socket.emit(mode === 'reasoning' ? 'reasoning_progress' : 'research_progress', {
                      type: mode === 'reasoning' ? 'reasoning_progress' : 'research_progress',
                      update: updateMsg,
                      data: { reasoning, sources, citations }
                    });
                  }
                }

                // Handle sources and citations
                if (parsedEvent.choices && parsedEvent.choices[0] && parsedEvent.choices[0].message) {
                    const message = parsedEvent.choices[0].message;
                    if (message.sources) {
                        sources = message.sources.map(s => ({ title: s.title || s.name, url: s.url, description: s.snippet }));
                         if (socket && mode === 'deep-research') {
                          socket.emit('research_progress', { 
                            type: 'research_progress', 
                            update: `Found ${sources.length} research sources.`, 
                            data: { sources, reasoning, citations } 
                          });
                        }
                    }
                    if (message.citations) {
                        citations = message.citations.map(c => ({ title: c.title || c.text, url: c.url }));
                         if (socket && mode === 'deep-research') {
                          socket.emit('research_progress', { 
                            type: 'research_progress', 
                            update: `Found ${citations.length} citations.`, 
                            data: { sources, reasoning, citations } 
                          });
                         }
                    }
                }

              } catch (e) {
                console.error('Error parsing JSON from SSE data:', e);
                console.error('Problematic JSON data string:', jsonData);
              }
            }
          }
        }
      }
      
      // Extract any remaining thinking process
      if (thinkingProcess === "") {
        thinkingProcess = extractThinkingProcess(fullContent);
      }
      
      // Filter out thinking tags from the final content
      fullContent = filterThinkingTags(fullContent);
      
      if (socket) {
        const finalUpdate = mode === 'reasoning' ? 'Analysis complete.' : 'Research complete.';
        socket.emit(mode === 'reasoning' ? 'reasoning_progress' : 'research_progress', {
          type: mode === 'reasoning' ? 'reasoning_progress' : 'research_progress',
          update: finalUpdate,
          data: { sources, citations, reasoning }
        });
      }
      
      return {
        content: fullContent || 'No streaming content received.',
        model: model,
        usage: { total_tokens: 0 },
        sources: sources,
        citations: citations,
        thinking: thinkingProcess,
        reasoning: reasoning
      };
    } else {
      throw new Error("Response body not available for streaming or not a streaming request.");
    }
  } catch (error) {
    console.error('Error in queryPerplexityStreaming:', error);
    if (error.name === 'AbortError') {
        throw new Error('Perplexity API request timed out.');
    }
    throw error;
  }
}

// Extract thinking process from content
function extractThinkingProcess(content) {
  if (!content) return "";
  
  // Extract <think>...</think> blocks
  const thinkMatch = /<think>([\s\S]*?)<\/think>/i.exec(content);
  if (thinkMatch && thinkMatch[1]) {
    return thinkMatch[1].trim();
  }
  
  // Extract <research>...</research> blocks
  const researchMatch = /<research>([\s\S]*?)<\/research>/i.exec(content);
  if (researchMatch && researchMatch[1]) {
    return researchMatch[1].trim();
  }
  
  return "";
}

// Function to filter out thinking tags and internal reasoning
function filterThinkingTags(content) {
  if (!content) return content;
  
  // Remove <think>...</think> blocks completely
  content = content.replace(/<think>[\s\S]*?<\/think>/gi, '');
  
  // Remove <research>...</research> blocks
  content = content.replace(/<research>[\s\S]*?<\/research>/gi, '');
  
  // Remove any remaining thinking artifacts
  content = content.replace(/^\s*thinking:[\s\S]*?(?=\n\n|\n[A-Z]|$)/im, '');
  content = content.replace(/^\s*analysis:[\s\S]*?(?=\n\n|\n[A-Z]|$)/im, '');
  content = content.replace(/^\s*research:[\s\S]*?(?=\n\n|\n[A-Z]|$)/im, '');
  
  // Clean up extra whitespace
  content = content.replace(/\n\s*\n\s*\n/g, '\n\n');
  content = content.trim();
  
  return content;
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
        const greetingMessage = "Hi! I'm Monday, your learning assistant powered by Perplexity Sonar. What would you like to dive into today?";
        
        console.log("Sending activation greeting:", greetingMessage);
        
        socket.emit('voice_response', {
          type: 'voice_response',
          message: greetingMessage,
          data: {
            panels: [],
            mode: 'basic',
            model: 'sonar',
            query: 'activation',
            citations: [],
            reasoning: [],
            sources: [],
            metadata: {
              tokensUsed: 0,
              responseTime: Date.now(),
              isActivation: true,
              youtubeVideoId: undefined
            }
          }
        });
        return;
      }

      // For all other cases (including Hey Monday + command), process with Perplexity
      const modeDetection = detectVoiceMode(command);
      const mode = modeDetection.mode;
      
      console.log(`Processing command in ${mode} mode (confidence: ${modeDetection.confidence}): "${command}"`);
      
      // For basic mode, send response immediately without any progress
      if (mode === 'basic') {
        try {
          console.log('Starting basic mode processing...');
          const perplexityResponse = await queryPerplexityStreaming(command, mode, null); // No socket for basic mode
          
          console.log(`Basic mode response ready: "${perplexityResponse.content.substring(0, 100)}..."`);
          console.log('Full response object:', {
            hasContent: !!perplexityResponse.content,
            contentLength: perplexityResponse.content?.length,
            model: perplexityResponse.model,
            usage: perplexityResponse.usage
          });

          // Create visualization data
          const visualizationData = {
            query: command,
            mode: mode,
            model: perplexityResponse.model,
            content: perplexityResponse.content,
            citations: perplexityResponse.citations || []
          };
          
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
          
          // Send response immediately for basic mode
          const responsePayload = {
            type: 'voice_response',
            message: perplexityResponse.content,
            data: {
              panels: panels,
              mode: mode,
              model: perplexityResponse.model,
              query: command,
              citations: [],
              reasoning: [],
              sources: [],
              visualizationData: visualizationData,
              metadata: {
                tokensUsed: perplexityResponse.usage?.total_tokens || 0,
                responseTime: Date.now()
              }
            }
          };
          
          console.log('Sending voice_response to client...');
          console.log('Response payload:', JSON.stringify(responsePayload).substring(0, 200) + '...');
          
          socket.emit('voice_response', responsePayload);
          
          console.log('Voice response sent successfully');
        } catch (error) {
          console.error('Error in basic mode:', error);
          console.error('Error stack:', error.stack);
          socket.emit('voice_error', {
            type: 'voice_error',
            error: `I encountered an issue: ${error.message}. Please try again.`
          });
        }
      } else {
        // For reasoning/research modes, use streaming with progress
        try {
          const perplexityResponse = await queryPerplexityStreaming(command, mode, socket);
          
          console.log(`${mode} mode response ready: "${perplexityResponse.content.substring(0, 100)}..."`);
          
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
            citations: perplexityResponse.citations || []
          }];
          
          // Send final response
          socket.emit('voice_response', {
            message: perplexityResponse.content,
            data: {
              panels: panels,
              mode: mode,
              model: perplexityResponse.model,
              query: command,
              citations: perplexityResponse.citations || [],
              reasoning: perplexityResponse.reasoning || [],
              sources: perplexityResponse.sources || [],
              thinking: perplexityResponse.thinking || '',
              metadata: {
                tokensUsed: perplexityResponse.usage?.total_tokens || 0,
                responseTime: Date.now()
              }
            }
          });
        } catch (error) {
          console.error(`Error in ${mode} mode:`, error);
          socket.emit('voice_error', {
            error: `I encountered an issue: ${error.message}. Please try again.`
          });
        }
      }
      
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