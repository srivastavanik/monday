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
const io = new Server(server, {
  cors: {
    origin: ["http://localhost:3000", "http://localhost:3001", "http://localhost:3002"],
    methods: ["GET", "POST"],
    credentials: true
  }
});

const PORT = 3001;

// Middleware
app.use(cors({
  origin: ["http://localhost:3000", "http://localhost:3001", "http://localhost:3002"],
  credentials: true
}));
app.use(express.json());

// Health endpoint
app.get('/health', (req, res) => {
  res.json({
    status: 'healthy',
    timestamp: new Date().toISOString(),
    uptime: process.uptime(),
    perplexityApiConfigured: !!process.env.PERPLEXITY_API_KEY
  });
});

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
            content: `You are Monday, an advanced AI learning companion and educational assistant. You have a friendly, knowledgeable personality and specialize in making complex topics easy to understand.

PERSONALITY:
- Enthusiastic about learning and teaching
- Patient and encouraging
- Uses clear, educational explanations
- Provides practical examples and analogies
- Encourages curiosity and deeper exploration

RESPONSE STYLE:
- Start responses naturally, as if speaking to a student
- Break down complex concepts into digestible parts
- Use analogies and real-world examples
- Suggest follow-up questions or related topics to explore
- Keep responses educational but engaging

CAPABILITIES:
- Explain any topic from basic to advanced levels
- Provide step-by-step tutorials
- Generate visual learning suggestions
- Recommend educational resources
- Answer follow-up questions with context

Always respond as Monday, the helpful learning companion. You are now in an active conversation with the user.`
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
  console.log('Client connected:', socket.id);

  socket.on('voice_command', async (data) => {
    console.log('Voice command received:', data);
    
    const { command, conversationActive, isExplicitTrigger, isActivation } = data;
    
    try {
      // If this is just the activation (Hey Monday without additional content), respond briefly
      if (isActivation && (!command || command.trim().length === 0)) {
        socket.emit('voice_response', {
          message: "Yes? I'm listening.",
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
      let mode = 'basic';
      const commandLower = command.toLowerCase();
      
      if (commandLower.includes('think about') || commandLower.includes('analyze')) {
        mode = 'reasoning';
      } else if (commandLower.includes('research into') || commandLower.includes('investigate')) {
        mode = 'deep-research';
      }
      
      console.log(`Processing command in ${mode} mode: "${command}"`);
      
      // Query Perplexity API with the actual command (removed progress update delay)
      const perplexityResponse = await queryPerplexity(command, mode);
      
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
          citations: [],
          reasoning: [],
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