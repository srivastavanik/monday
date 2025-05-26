import axios from 'axios'
import { logger } from '../utils/logger.js'

interface ConversationEntry {
  role: 'user' | 'assistant';
  content: string;
  ttsContent?: string;
  timestamp: number;
}

export interface PerplexityQuery {
  query: string
  mode: 'basic' | 'reasoning' | 'research'
  context?: string[]
  sessionId?: string
  progressCallback?: (update: ProgressUpdate) => void
}

export interface ProgressUpdate {
  type: 'thinking' | 'researching' | 'searching' | 'analyzing' | 'synthesizing' | 'complete'
  message: string
  progress: number // 0-100
  sources?: string[]
  reasoning?: string
  metadata?: any
}

export interface PerplexityResponse {
  id: string
  model: string
  content: string
  fullContent?: string
  finalPreview?: string
  citations?: Citation[]
  reasoning?: ReasoningStep[]
  sources?: Source[]
  metadata: {
    tokensUsed: number
    responseTime: number
    confidence?: number
    isThinking?: boolean
    isResearching?: boolean
    isStreaming?: boolean
    progressUpdates?: ProgressUpdate[]
  }
}

export interface Citation {
  id: string
  url: string
  title: string
  snippet: string
  publishedDate?: string
  domain: string
}

export interface ReasoningStep {
  step: number
  content: string
  confidence: number
  sources: string[]
  timestamp?: number
}

export interface Source {
  id: string
  url: string
  title: string
  snippet: string
  relevanceScore: number
  publishedDate?: string
}

interface QueryContext {
  service?: string;
  isMondayActivation?: boolean;
  isInActiveConversation?: boolean;
  sessionInConversation?: boolean;
}

// Define ContextCleaner class locally to avoid import issues
class ContextCleaner {
  // Remove voice recognition artifacts and clean messages
  static cleanMessage(message: string): string {
    // Remove common artifacts
    const artifacts = [
      /^s for you to explore$/i,
      /^okay$/i,
      /^mm-hmm$/i,
      /^uh$/i,
      /^\w{1,2}$/i, // Single or double character fragments
      /^thanks for sharing your curiosity/i, // Monday's own phrases being picked up
      /^i've gathered more details/i,
      /^that's a great topic/i,
    ];
    
    for (const artifact of artifacts) {
      if (artifact.test(message.trim())) {
        return ''; // Return empty to be filtered out
      }
    }
    
    // Clean up truncated sentences
    if (message.endsWith('...') || message.endsWith('…')) {
      // This is fine, keep it
    } else if (message.length < 10 && !message.endsWith('.') && !message.endsWith('?') && !message.endsWith('!')) {
      // Likely a fragment
      return '';
    }
    
    return message.trim();
  }
  
  static validateAndCleanContext(entries: ConversationEntry[]): Array<{role: 'user' | 'assistant', content: string}> {
    const cleaned: Array<{role: 'user' | 'assistant', content: string}> = [];
    
    for (const entry of entries) {
      const cleanContent = this.cleanMessage(entry.content);
      
      // Skip empty or invalid messages
      if (!cleanContent || cleanContent.length < 3) {
        continue;
      }
      
      cleaned.push({
        role: entry.role,
        content: cleanContent
      });
    }
    
    return cleaned;
  }
  
  // Convert old string format to new format for backward compatibility
  static convertLegacyContext(legacyContext: string[]): ConversationEntry[] {
    const entries: ConversationEntry[] = [];
    
    for (const ctx of legacyContext) {
      if (ctx.startsWith('User: ')) {
        entries.push({
          role: 'user',
          content: ctx.replace('User: ', ''),
          timestamp: Date.now()
        });
      } else if (ctx.startsWith('Monday: ')) {
        entries.push({
          role: 'assistant',
          content: ctx.replace('Monday: ', ''),
          timestamp: Date.now()
        });
      }
    }
    
    return entries;
  }
}

class PerplexityService {
  private apiKey: string
  private readonly baseUrl: string = 'https://api.perplexity.ai'
  private conversationHistory: string[] = []
  private lastRequestTime: number = 0
  private requestCount: number = 0
  private readonly MIN_REQUEST_INTERVAL = 2000 // 2 seconds between requests
  private readonly MAX_REQUESTS_PER_MINUTE = 20 // Limit to 20 requests per minute
  private requestTimes: number[] = []
  private systemPrompt = `You are Monday, an advanced AI learning companion for VR education, powered by Perplexity Sonar.

Core Identity:
- Intelligent, curious, and passionate about learning
- Speak conversationally and encouragingly 
- Keep responses clear and TTS-friendly (avoid excessive symbols)
- Show genuine interest in helping users learn

Your Role:
- Guide users through immersive learning experiences
- Provide clear, educational responses with context
- Encourage deeper exploration of topics
- Mention when topics might benefit from reasoning or research modes

Response Guidelines:
- Keep responses conversational (2-4 sentences for basic queries)
- Use natural speech patterns suitable for voice synthesis
- Always end with engagement (questions, suggestions, or offers to explore more)
- Keep responses concise and complete - avoid getting cut off`
  
  constructor() {
    this.apiKey = process.env.PERPLEXITY_API_KEY || ''
    if (!this.apiKey) {
      throw new Error('PERPLEXITY_API_KEY environment variable is not set')
    }
    console.log('[DEBUG] API Key length:', this.apiKey.length)
    console.log('[DEBUG] API Key first 10 chars:', this.apiKey.substring(0, 10))
    console.log('[DEBUG] Base URL:', this.baseUrl)
  }

  private async makeRequest(endpoint: string, data: any): Promise<any> {
    // RATE LIMITING: Check if we're making too many requests
    const now = Date.now()
    
    // Remove requests older than 1 minute
    this.requestTimes = this.requestTimes.filter(time => now - time < 60000)
    
    // Check if we've exceeded the rate limit
    if (this.requestTimes.length >= this.MAX_REQUESTS_PER_MINUTE) {
      const oldestRequest = Math.min(...this.requestTimes)
      const waitTime = 60000 - (now - oldestRequest)
      console.warn(`[RATE LIMIT] Too many requests. Waiting ${waitTime}ms before next request.`)
      await new Promise(resolve => setTimeout(resolve, waitTime))
    }
    
    // Check minimum interval between requests
    const timeSinceLastRequest = now - this.lastRequestTime
    if (timeSinceLastRequest < this.MIN_REQUEST_INTERVAL) {
      const waitTime = this.MIN_REQUEST_INTERVAL - timeSinceLastRequest
      console.log(`[RATE LIMIT] Waiting ${waitTime}ms to maintain minimum interval.`)
      await new Promise(resolve => setTimeout(resolve, waitTime))
    }
    
    // Record this request
    this.requestTimes.push(Date.now())
    this.lastRequestTime = Date.now()
    this.requestCount++
    
    console.log(`[RATE LIMIT] Request #${this.requestCount}, ${this.requestTimes.length} requests in last minute`)
    
    const startTime = Date.now()
    
    // Ensure endpoint starts with /
    const cleanEndpoint = endpoint.startsWith('/') ? endpoint : `/${endpoint}`
    const fullUrl = `${this.baseUrl}${cleanEndpoint}`
    
    try {
      // Log the request details for debugging
      console.log('[DEBUG] Making request to Perplexity API:', {
        baseUrl: this.baseUrl,
        endpoint: cleanEndpoint,
        fullUrl: fullUrl,
        apiKey: this.apiKey.substring(0, 10) + '...',
        dataKeys: Object.keys(data)
      })

      const response = await axios({
        method: 'post',
        url: fullUrl,
        data: data,
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${this.apiKey.trim()}`,
          'Accept': 'application/json'
        },
        timeout: 60000 // Increased timeout for streaming
      })
      
      const responseTime = Date.now() - startTime
      
      // Log the complete response for debugging
      console.log('[DEBUG] Perplexity API response:', {
        status: response.status,
        hasData: !!response.data,
        dataKeys: response.data ? Object.keys(response.data) : []
      })
      
      logger.info('Perplexity API request completed', {
        endpoint: cleanEndpoint,
        responseTime: `${responseTime}ms`,
        statusCode: response.status
      })
      
      return response.data
    } catch (error: any) {
      const responseTime = Date.now() - startTime
      console.error('[DEBUG] Perplexity API request failed:', {
        fullUrl: fullUrl,
        error: error.message,
        status: error.response?.status,
        hostname: error.hostname,
        code: error.code
      })
      
      logger.error('Perplexity API request failed', {
        endpoint: cleanEndpoint,
        responseTime: `${responseTime}ms`,
        error: error.message,
        status: error.response?.status,
        responseData: error.response?.data ? JSON.stringify(error.response.data) : 'No response data'
      })
      throw error
    }
  }

  private async makeStreamingRequest(endpoint: string, data: any, progressCallback?: (update: ProgressUpdate) => void): Promise<any> {
    const startTime = Date.now()
    const cleanEndpoint = endpoint.startsWith('/') ? endpoint : `/${endpoint}`
    const fullUrl = `${this.baseUrl}${cleanEndpoint}`
    
    try {
      console.log('[DEBUG] Making streaming request to Perplexity API:', {
        fullUrl: fullUrl,
        model: data.model,
        streaming: data.stream,
        maxTokens: data.max_tokens
      })

      const response = await fetch(fullUrl, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${this.apiKey.trim()}`,
          'Accept': 'text/event-stream'
        },
        body: JSON.stringify(data)
      })

      if (!response.ok) {
        const errorText = await response.text()
        console.error('[DEBUG] Streaming request failed:', {
          status: response.status,
          statusText: response.statusText,
          errorText: errorText
        })
        throw new Error(`HTTP ${response.status}: ${response.statusText} - ${errorText}`)
      }

      if (!response.body) {
        throw new Error('No response body for streaming')
      }

      const reader = response.body.getReader()
      const decoder = new TextDecoder()
      let fullContent = ''
      let buffer = ''
      let progressCount = 0
      let lastProgressUpdate = 0
      let isComplete = false

      try {
        while (true) {
          const { done, value } = await reader.read()
          
          if (done) {
            console.log('[DEBUG] Streaming completed, final content length:', fullContent.length)
            break
          }

          buffer += decoder.decode(value, { stream: true })
          const lines = buffer.split('\n')
          buffer = lines.pop() || '' // Keep incomplete line in buffer

          for (const line of lines) {
            if (line.trim() === '') continue
            if (line.startsWith('data: ')) {
              const sseData = line.slice(6)
              if (sseData === '[DONE]') {
                console.log('[DEBUG] Received [DONE] signal, content length:', fullContent.length)
                isComplete = true
                continue
              }

              try {
                const parsed = JSON.parse(sseData)
                const delta = parsed.choices?.[0]?.delta?.content
                
                if (delta) {
                  fullContent += delta
                  progressCount++

                  // IMPROVED progress updates - send every 3 chunks AND at least 800ms apart
                  const now = Date.now()
                  if (progressCallback && progressCount % 3 === 0 && (now - lastProgressUpdate) >= 800) {
                    const progress = Math.min(85, Math.floor(progressCount * 2)) // Moderate progress, cap at 85%
                    
                    let updateType: ProgressUpdate['type'] = 'thinking'
                    let message = 'Processing your request...'
                    
                    if (data.model?.includes('reasoning')) {
                      updateType = 'thinking'
                      message = 'Working through the reasoning process...'
                    } else if (data.model?.includes('research')) {
                      updateType = 'researching'
                      message = 'Gathering information from multiple sources...'
                    }

                    console.log(`[DEBUG] Sending progress update: ${progress}% (chunk ${progressCount}, content: ${fullContent.length} chars)`)
                    progressCallback({
                      type: updateType,
                      message: message,
                      progress: progress,
                      reasoning: fullContent // Send full content for real-time display
                    })
                    
                    lastProgressUpdate = now
                  }
                }
              } catch (parseError) {
                console.warn('Failed to parse streaming chunk:', parseError)
              }
            }
          }
        }
      } finally {
        reader.releaseLock()
      }

      // Ensure we have complete content
      console.log('[DEBUG] Final content check:', {
        contentLength: fullContent.length,
        isComplete: isComplete,
        lastChars: fullContent.slice(-50)
      })

      // Send completion update with full content
      if (progressCallback) {
        console.log('[DEBUG] Sending completion update with full content')
        progressCallback({
          type: 'complete',
          message: 'Analysis complete!',
          progress: 100,
          reasoning: fullContent
        })
      }

      const responseTime = Date.now() - startTime
      console.log('[DEBUG] Streaming request completed:', {
        responseTime: `${responseTime}ms`,
        contentLength: fullContent.length,
        totalChunks: progressCount,
        isComplete: isComplete
      })

      // Return in expected format
      return {
        id: `streaming_${Date.now()}`,
        model: data.model,
        choices: [{
          message: {
            content: fullContent
          }
        }],
        usage: {
          total_tokens: Math.ceil(fullContent.length / 4) // Rough estimate
        }
      }

    } catch (error: any) {
      const responseTime = Date.now() - startTime
      console.error('[DEBUG] Streaming request failed:', {
        fullUrl: fullUrl,
        error: error.message,
        responseTime: `${responseTime}ms`
      })
      throw error
    }
  }

  public async basicQuery(query: string, context: QueryContext): Promise<PerplexityResponse> {
    try {
      // Format conversation history with alternating roles
      const formattedHistory = this.conversationHistory.map((msg, index) => ({
        role: index % 2 === 0 ? 'user' : 'assistant',
        content: msg
      }))

      // Add current query
      const messages = [
        {
          role: 'system',
          content: this.systemPrompt
        },
        ...formattedHistory,
        {
          role: 'user',
          content: query
        }
      ]

      const response = await this.makeRequest('/chat/completions', {
        model: 'sonar',
        messages: messages,
        max_tokens: 150, // REDUCED from 300 to prevent excessive costs
        temperature: 0.7,
        top_p: 0.9,
        stream: false
      })

      // Create short TTS response and full content
      const fullContent = response.choices[0].message.content
      const shortResponse = this.createShortTTSResponse(fullContent, query)

      return {
        id: response.id,
        model: response.model,
        content: shortResponse,
        fullContent: fullContent,
        citations: this.extractCitations(response),
        metadata: {
          tokensUsed: response.usage?.total_tokens || 0,
          responseTime: Date.now() - Date.now()
        }
      }
    } catch (error) {
      console.error('Perplexity basicQuery error:', error)
      throw error
    }
  }

  public async reasoningQuery(query: string, context?: string[] | ConversationEntry[], progressCallback?: (update: ProgressUpdate) => void): Promise<PerplexityResponse> {
    console.log('🔥 NEW STREAMING REASONING QUERY METHOD CALLED!', { query, contextLength: context?.length });
    
    const messages: any[] = [];
    const progressUpdates: ProgressUpdate[] = [];
    
    // Add system message
    const systemPrompt = `You are Monday, an AI learning companion with advanced reasoning capabilities.

Your Reasoning Approach:
- Break down complex problems into clear, logical steps
- Show your thinking process transparently
- Provide confidence levels for each reasoning step
- Connect concepts and show relationships
- Use analogies and examples to make concepts accessible

Response Format:
- Start with a brief overview of your approach
- Present 3-5 clear reasoning steps
- Each step should be conversational and TTS-friendly
- End with synthesis and suggestions for further exploration
- IMPORTANT: Always complete your thoughts - don't leave sentences unfinished
- Structure your response to fit within the token limit while being complete

Educational Focus:
- Help users understand not just what, but why and how
- Encourage critical thinking
- Make complex topics approachable
- Ensure every response has a clear conclusion`;

    messages.push({
      role: 'system',
      content: systemPrompt
    });
    
    console.log('🔥 Processing context...', { hasContext: !!context, contextType: typeof context?.[0] });
    
    // Process context with new structure
    if (context && Array.isArray(context)) {
      let structuredContext: ConversationEntry[];
      
      // Handle both old string format and new ConversationEntry format
      if (context.length > 0 && typeof context[0] === 'string') {
        console.log('🔥 Converting legacy context format');
        structuredContext = ContextCleaner.convertLegacyContext(context as string[]);
      } else {
        console.log('🔥 Using new context format');
        structuredContext = context as ConversationEntry[];
      }
      
      console.log('🔥 Structured context:', structuredContext);
      
      // Clean the context
      const cleanedMessages = ContextCleaner.validateAndCleanContext(structuredContext);
      console.log('🔥 Cleaned messages:', cleanedMessages);
      
      // Ensure perfect alternation
      const alternatingMessages = this.ensureAlternation(cleanedMessages);
      console.log('🔥 Alternating messages:', alternatingMessages);
      
      // Add to messages
      messages.push(...alternatingMessages);
    }
    
    // Add current query
    messages.push({
      role: 'user',
      content: `Please think through this step by step: ${query}`
    });
    
    console.log('🔥 Final messages before validation:', messages);
    
    // Final validation
    this.validateMessageStructure(messages);
    
    // Log for debugging
    this.logApiRequest('/chat/completions', {
      model: 'sonar-reasoning-pro',
      messages: messages
    });
    
    const requestData = {
      model: 'sonar-reasoning-pro',
      messages: messages,
      max_tokens: 600, // INCREASED from 400 to ensure complete responses
      temperature: 0.2,
      stream: true // Enable streaming for progressive display
    };

    // Send initial progress update
    if (progressCallback) {
      const initialUpdate: ProgressUpdate = {
        type: 'thinking',
        message: `I'm thinking through ${query.replace(/^(please\s+)?think\s+(through\s+|about\s+)?/i, '').trim()} step by step. Let me work through this systematically.`,
        progress: 10
      };
      progressUpdates.push(initialUpdate);
      progressCallback(initialUpdate);
    }

    try {
      const result = await this.makeStreamingRequest('/chat/completions', requestData, (update) => {
        // Override update type for reasoning
        const reasoningUpdate: ProgressUpdate = {
          ...update,
          type: update.type === 'thinking' ? 'synthesizing' : update.type,
          message: update.type === 'thinking' ? 'Synthesizing reasoning findings...' : update.message
        };
        progressUpdates.push(reasoningUpdate);
        if (progressCallback) progressCallback(reasoningUpdate);
      });
      
      const fullContent = result.choices?.[0]?.message?.content || 'No response generated';
      
      // Check if content appears to be truncated
      const seemsTruncated = fullContent.length > 50 && (
        !fullContent.trim().endsWith('.') && 
        !fullContent.trim().endsWith('!') && 
        !fullContent.trim().endsWith('?') &&
        !fullContent.trim().endsWith('</think>')
      );
      
      if (seemsTruncated) {
        console.warn('[DEBUG] ⚠️ REASONING CONTENT MAY BE TRUNCATED:', {
          contentLength: fullContent.length,
          lastChars: fullContent.slice(-100),
          tokensUsed: result.usage?.total_tokens || 0,
          maxTokens: requestData.max_tokens
        });
      } else {
        console.log('[DEBUG] ✅ Reasoning content appears complete:', {
          contentLength: fullContent.length,
          tokensUsed: result.usage?.total_tokens || 0,
          maxTokens: requestData.max_tokens
        });
      }
      
      // Create a thinking message for TTS and main panel
      const topicName = query.replace(/^(please\s+)?think\s+(through\s+|about\s+)?/i, '').trim();
      const thinkingMessage = `I'm going to think through ${topicName} step by step. Let me break this down systematically and analyze the different aspects, approaches, and key considerations. I'll work through this methodically to give you a comprehensive understanding.`;
      
      console.log('🎯 Generated thinking message:', {
        topicName: topicName,
        messageLength: thinkingMessage.length,
        message: thinkingMessage.substring(0, 100)
      });
      
      // Generate a final response preview from the reasoning content
      let finalPreview = '';
      if (fullContent && fullContent.length > 100) {
        // Extract key insights from the reasoning for the preview
        const cleanReasoning = fullContent.replace(/<think>|<\/think>/g, '').trim();
        const sentences = cleanReasoning.split(/[.!?]+/).filter(s => s.trim().length > 20);
        
        if (sentences.length >= 3) {
          // Take key sentences and create a summary
          const keyPoints = sentences.slice(0, 3).map(s => s.trim()).join('. ');
          finalPreview = `Based on my analysis: ${keyPoints}. I've broken this down into clear steps for you to explore.`;
        } else {
          finalPreview = `I've analyzed ${query.replace(/^(please\s+)?think\s+(through\s+|about\s+)?/i, '').trim()} and broken it down into logical steps. The reasoning process reveals key insights about the different approaches and their applications.`;
        }
      } else {
        finalPreview = `I've completed my analysis of ${query.replace(/^(please\s+)?think\s+(through\s+|about\s+)?/i, '').trim()}. The reasoning process is now available for you to review.`;
      }
      
      console.log('🎯 Generated final preview:', {
        previewLength: finalPreview.length,
        preview: finalPreview.substring(0, 100)
      });
      
      return {
        id: result.id || 'reasoning_query',
        model: result.model || 'sonar-reasoning-pro',
        content: thinkingMessage, // Short thinking message for TTS
        fullContent: fullContent, // Full reasoning for progressive display
        finalPreview: finalPreview, // Final response preview for TTS after completion
        citations: this.extractCitations(result),
        reasoning: this.extractReasoningSteps(fullContent),
        metadata: {
          tokensUsed: result.usage?.total_tokens || 0,
          responseTime: 0,
          isThinking: true, // Flag to indicate this is a thinking response
          isStreaming: true,
          progressUpdates: progressUpdates
        }
      };
    } catch (error: any) {
      console.error('Perplexity API Error:', error.response?.data || error.message);
      throw error;
    }
  }

  public async deepResearch(query: string, context?: string[] | ConversationEntry[], progressCallback?: (update: ProgressUpdate) => void): Promise<PerplexityResponse> {
    console.log('🔥 NEW STREAMING DEEP RESEARCH METHOD CALLED!', { query, contextLength: context?.length });
    
    const messages: any[] = [];
    const progressUpdates: ProgressUpdate[] = [];
    
    // Add system message
    const systemPrompt = `You are Monday, conducting comprehensive research analysis.

Research Methodology:
- Synthesize information from multiple high-quality sources
- Present multiple perspectives on complex topics
- Evaluate source credibility and recency
- Identify knowledge gaps and areas of debate
- Connect findings to broader implications

Response Structure:
- Opening: Brief context and research scope
- Main findings: 3-4 key insights with source backing
- Analysis: Critical evaluation and synthesis
- Implications: Broader significance and applications
- Conclusion: Summary and further research directions
- IMPORTANT: Always complete your analysis - don't leave thoughts unfinished

Voice-Friendly Delivery:
- Use clear, flowing language suitable for TTS
- Break up long sections with natural pauses
- Avoid excessive technical jargon without explanation
- Maintain conversational tone despite depth
- Ensure every response has a clear conclusion`;

    messages.push({
      role: 'system',
      content: systemPrompt
    });

    console.log('🔥 Processing research context...', { hasContext: !!context, contextType: typeof context?.[0] });

    // Process context with new structure
    if (context && Array.isArray(context)) {
      let structuredContext: ConversationEntry[];
      
      // Handle both old string format and new ConversationEntry format
      if (context.length > 0 && typeof context[0] === 'string') {
        console.log('🔥 Converting legacy research context format');
        structuredContext = ContextCleaner.convertLegacyContext(context as string[]);
      } else {
        console.log('🔥 Using new research context format');
        structuredContext = context as ConversationEntry[];
      }
      
      console.log('🔥 Structured research context:', structuredContext);
      
      // Clean the context
      const cleanedMessages = ContextCleaner.validateAndCleanContext(structuredContext);
      console.log('🔥 Cleaned research messages:', cleanedMessages);
      
      // Ensure perfect alternation
      const alternatingMessages = this.ensureAlternation(cleanedMessages);
      console.log('🔥 Alternating research messages:', alternatingMessages);
      
      // Add to messages
      messages.push(...alternatingMessages);
    }

    // Add the user query
    messages.push({
      role: 'user',
      content: `Conduct a comprehensive research analysis on: ${query}`
    });

    console.log('🔥 Final research messages before validation:', messages);

    // Final validation
    this.validateMessageStructure(messages);
    
    // Log for debugging
    this.logApiRequest('/chat/completions', {
      model: 'sonar-deep-research',
      messages: messages
    });

    // Send initial progress update
    if (progressCallback) {
      const initialUpdate: ProgressUpdate = {
        type: 'researching',
        message: `I'm conducting comprehensive research on ${query.replace(/^(please\s+)?(research\s+|investigate\s+)?/i, '').trim()}. Let me gather information from multiple sources and analyze the findings.`,
        progress: 5,
        sources: []
      };
      progressUpdates.push(initialUpdate);
      progressCallback(initialUpdate);

      // Simulate research phases with progress updates
      setTimeout(() => {
        const searchUpdate: ProgressUpdate = {
          type: 'searching',
          message: 'Searching through academic databases and reliable sources...',
          progress: 25,
          sources: ['Academic databases', 'Research papers', 'Expert publications']
        };
        progressUpdates.push(searchUpdate);
        progressCallback(searchUpdate);
      }, 1000);

      setTimeout(() => {
        const analyzeUpdate: ProgressUpdate = {
          type: 'analyzing',
          message: 'Analyzing source credibility and synthesizing findings...',
          progress: 60,
          sources: ['Peer-reviewed sources', 'Recent publications', 'Expert opinions']
        };
        progressUpdates.push(analyzeUpdate);
        progressCallback(analyzeUpdate);
      }, 3000);
    }

    const requestData = {
      model: 'sonar-deep-research',
      messages: messages,
      max_tokens: 700, // INCREASED from 500 to ensure complete responses
      temperature: 0.3,
      stream: true // Enable streaming for progressive display
    };

    try {
      const result = await this.makeStreamingRequest('/chat/completions', requestData, (update) => {
        // Override update type for research
        const researchUpdate: ProgressUpdate = {
          ...update,
          type: update.type === 'thinking' ? 'synthesizing' : update.type,
          message: update.type === 'thinking' ? 'Synthesizing research findings...' : update.message
        };
        progressUpdates.push(researchUpdate);
        if (progressCallback) progressCallback(researchUpdate);
      });
      
      const fullContent = result.choices?.[0]?.message?.content || 'No response generated';
      
      // Check if content appears to be truncated
      const seemsTruncated = fullContent.length > 50 && (
        !fullContent.trim().endsWith('.') && 
        !fullContent.trim().endsWith('!') && 
        !fullContent.trim().endsWith('?')
      );
      
      if (seemsTruncated) {
        console.warn('[DEBUG] ⚠️ RESEARCH CONTENT MAY BE TRUNCATED:', {
          contentLength: fullContent.length,
          lastChars: fullContent.slice(-100),
          tokensUsed: result.usage?.total_tokens || 0,
          maxTokens: requestData.max_tokens
        });
      } else {
        console.log('[DEBUG] ✅ Research content appears complete:', {
          contentLength: fullContent.length,
          tokensUsed: result.usage?.total_tokens || 0,
          maxTokens: requestData.max_tokens
        });
      }
      
      // Create a research message for TTS and main panel
      const topicName = query.replace(/^(please\s+)?(research\s+|investigate\s+)?/i, '').trim();
      const researchMessage = `I'm going to conduct comprehensive research on ${topicName}. Let me gather information from multiple high-quality sources, analyze different perspectives, and synthesize the findings. I'll examine the latest developments, key insights, and provide you with a thorough analysis.`;
      
      console.log('🎯 Generated research message:', {
        topicName: topicName,
        messageLength: researchMessage.length,
        message: researchMessage.substring(0, 100)
      });
      
      return {
        id: result.id || 'research_query',
        model: result.model || 'sonar-deep-research',
        content: researchMessage, // Short research message for TTS
        fullContent: fullContent, // Full research for progressive display
        citations: this.extractCitations(result),
        sources: this.extractSources(result),
        metadata: {
          tokensUsed: result.usage?.total_tokens || 0,
          responseTime: 0,
          isResearching: true, // Flag to indicate this is a research response
          isStreaming: true,
          progressUpdates: progressUpdates
        }
      };
    } catch (error: any) {
      console.error('Perplexity API Error:', error.response?.data || error.message);
      throw error;
    }
  }

  private ensureAlternation(messages: Array<{role: string, content: string}>): Array<{role: string, content: string}> {
    if (messages.length === 0) return [];
    
    const result: Array<{role: string, content: string}> = [];
    let expectedRole: 'user' | 'assistant' = 'user'; // First message after system should be user
    
    for (const msg of messages) {
      if (msg.role === 'system') {
        continue; // Skip system messages in this logic
      }
      
      if (msg.role === expectedRole) {
        result.push(msg);
        // Toggle expected role
        expectedRole = expectedRole === 'user' ? 'assistant' : 'user';
      } else {
        // Wrong role order detected - skip this message but log it
        console.warn(`⚠️ Skipping message with role ${msg.role}, expected ${expectedRole}: "${msg.content.substring(0, 50)}..."`);
      }
    }
    
    // CRITICAL FIX: Ensure we end with assistant message so the new user query creates proper alternation
    // If we have messages and the last one is user, we need to add a placeholder assistant message
    if (result.length > 0 && result[result.length - 1].role === 'user') {
      result.push({
        role: 'assistant',
        content: 'I understand. Please continue.'
      });
      console.log('🔧 Added placeholder assistant message to ensure alternation');
    }
    
    console.log('🔧 Final alternation result:', result.map(m => `${m.role}: ${m.content.substring(0, 30)}...`));
    
    return result;
  }

  private validateMessageStructure(messages: any[]): void {
    let lastRole: string | null = null;
    
    for (let i = 0; i < messages.length; i++) {
      const msg = messages[i];
      
      // System messages can be at the start
      if (msg.role === 'system' && i > 0 && messages[i-1].role !== 'system') {
        throw new Error('System messages must be at the beginning');
      }
      
      // After system messages, roles must alternate
      if (msg.role !== 'system') {
        if (lastRole && lastRole !== 'system' && msg.role === lastRole) {
          throw new Error(`Role alternation violated at index ${i}: ${lastRole} -> ${msg.role}`);
        }
        lastRole = msg.role;
      }
    }
    
    // Must end with user message (the current query)
    if (messages[messages.length - 1].role !== 'user') {
      throw new Error('Messages must end with user role');
    }
  }

  private logApiRequest(endpoint: string, requestBody: any): void {
    console.log('\n=== PERPLEXITY API REQUEST ===');
    console.log('Endpoint:', endpoint);
    console.log('Model:', requestBody.model);
    console.log('Message count:', requestBody.messages.length);
    console.log('Streaming:', requestBody.stream || false);
    console.log('Messages:');
    
    requestBody.messages.forEach((msg: any, index: number) => {
      console.log(`[${index}] Role: ${msg.role}, Content: ${msg.content.substring(0, 50)}...`);
    });
    
    // Validate alternation
    let lastRole = null;
    let alternationValid = true;
    
    for (let i = 0; i < requestBody.messages.length; i++) {
      const msg = requestBody.messages[i];
      if (msg.role !== 'system') {
        if (lastRole && lastRole === msg.role) {
          console.error(`❌ ALTERNATION ERROR at index ${i}: ${lastRole} -> ${msg.role}`);
          alternationValid = false;
        }
        lastRole = msg.role;
      }
    }
    
    console.log('Alternation valid:', alternationValid ? '✅' : '❌');
    console.log('============================\n');
  }

  private extractCitations(response: any): Citation[] {
    if (!response.citations) return []
    
    return response.citations.map((citation: any, index: number) => ({
      id: `citation_${index}`,
      url: citation.url || citation,
      title: citation.title || `Source ${index + 1}`,
      snippet: citation.snippet || citation.text || ''
    }))
  }

  private extractReasoningSteps(content: string): ReasoningStep[] {
    const steps: ReasoningStep[] = []
    const lines = content.split('\n')
    let stepCount = 0
    
    for (const line of lines) {
      // Look for step indicators like "Step 1:", "1.", etc.
      const stepMatch = line.match(/^(?:Step\s+)?(\d+)[:.]?\s*(.+)$/i)
      if (stepMatch) {
        stepCount++
        steps.push({
          step: stepCount,
          content: stepMatch[2].trim(),
          confidence: 0.8, // Default confidence
          sources: [],
          timestamp: Date.now()
        })
      }
    }
    
    return steps
  }

  private extractSources(result: any): Source[] {
    const sources: Source[] = []
    
    if (result.citations) {
      result.citations.forEach((citation: any, index: number) => {
        sources.push({
          id: `source_${index}`,
          url: citation.url || '',
          title: citation.title || 'Untitled',
          snippet: citation.text || citation.snippet || '',
          relevanceScore: 0.8, // Default relevance
          publishedDate: citation.published_date
        })
      })
    }
    
    return sources
  }

  private extractDomain(url: string): string {
    try {
      return new URL(url).hostname
    } catch {
      return 'unknown'
    }
  }

  private createShortTTSResponse(fullContent: string, query: string): string {
    console.log('🔥 TTS PROCESSING:', {
      fullContentPreview: fullContent?.substring(0, 100),
      fullContentLength: fullContent?.length,
      query: query
    });
    
    // Create a proper introductory sentence for TTS that doesn't cut off
    const sentences = fullContent.split(/[.!?]+/).filter(s => s.trim().length > 0)
    
    console.log('🔥 TTS SENTENCES:', {
      sentenceCount: sentences.length,
      firstSentence: sentences[0]?.substring(0, 100)
    });
    
    if (sentences.length === 0) {
      const fallback = "I found some information about that topic for you.";
      console.log('🔥 TTS FALLBACK (no sentences):', fallback);
      return fallback;
    }
    
    // Take the first complete sentence and ensure it's a good TTS intro
    let firstSentence = sentences[0].trim()
    
    // If the first sentence is too long (>120 chars), create a custom intro
    if (firstSentence.length > 120) {
      // Extract key topic from the query for a personalized intro
      const cleanQuery = query.replace(/^(hey monday,?\s*)/i, '').trim()
      if (cleanQuery.length > 0) {
        const customResponse = `I found some great information about ${cleanQuery}. Let me share what I discovered.`;
        console.log('🔥 TTS CUSTOM (long sentence):', customResponse);
        return customResponse;
      } else {
        const fallback = "I found some interesting information to share with you.";
        console.log('🔥 TTS FALLBACK (long sentence, no query):', fallback);
        return fallback;
      }
    }
    
    // Ensure the sentence ends properly
    if (!firstSentence.endsWith('.') && !firstSentence.endsWith('!') && !firstSentence.endsWith('?')) {
      firstSentence += '.'
    }
    
    // Add a connecting phrase to indicate there's more detail in the panel
    const finalResponse = `${firstSentence} I've gathered more details for you to explore.`;
    console.log('🔥 TTS FINAL RESPONSE:', finalResponse);
    return finalResponse;
  }

  async processQuery(queryData: PerplexityQuery): Promise<PerplexityResponse> {
    try {
      const queryContext: QueryContext = {
        service: 'monday-backend',
        isMondayActivation: queryData.query.toLowerCase().includes('hey monday'),
        isInActiveConversation: queryData.context && queryData.context.length > 0,
        sessionInConversation: queryData.context && queryData.context.length > 0
      }
      
      return await this.basicQuery(queryData.query, queryContext)
    } catch (error) {
      console.error('[DEBUG] Failed to process Perplexity query:', error)
      throw error
    }
  }
}

export const perplexityService = new PerplexityService() 