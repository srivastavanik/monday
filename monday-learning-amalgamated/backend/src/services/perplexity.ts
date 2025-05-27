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
  thinkingProcess?: string
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
    ttsMessage?: string
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

  private async makeStreamingRequest(endpoint: string, data: any, progressCallback?: (update: any) => void): Promise<any> {
    console.log('[DEBUG] Making streaming request to Perplexity API:', {
      fullUrl: `${this.baseUrl}${endpoint}`,
      model: data.model,
      streaming: data.stream,
      maxTokens: data.max_tokens
    });

    const response = await fetch(`${this.baseUrl}${endpoint}`, {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${this.apiKey}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify(data),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`HTTP ${response.status}: ${errorText}`);
    }

    if (!response.body) {
      throw new Error('No response body available for streaming');
    }

    const reader = response.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    let fullContent = '';
    let chunkCount = 0;
    let lastProgressUpdate = 0;
    const startTime = Date.now();

    try {
      while (true) {
        const { done, value } = await reader.read();
        
        if (done) {
          console.log('[DEBUG] Streaming completed, final content length:', fullContent.length);
          break;
        }

        buffer += decoder.decode(value, { stream: true });
        const lines = buffer.split('\n');
        buffer = lines.pop() || '';

        for (const line of lines) {
          if (line.startsWith('data: ')) {
            const jsonStr = line.slice(6);
            if (jsonStr === '[DONE]') {
              console.log('[DEBUG] Received [DONE] signal');
              continue;
            }

            try {
              const chunk = JSON.parse(jsonStr);
              const content = chunk.choices?.[0]?.delta?.content || '';
              
              if (content) {
                fullContent += content;
                chunkCount++;
                
                // Send progress updates every 10 chunks or every 2 seconds, whichever comes first
                const now = Date.now();
                const shouldSendUpdate = (
                  chunkCount % 10 === 0 || 
                  (now - lastProgressUpdate) >= 2000 ||
                  fullContent.length > 0 && fullContent.length % 500 === 0
                );
                
                if (shouldSendUpdate && progressCallback) {
                  // Calculate progress based on content length and estimated completion
                  let estimatedProgress = Math.min(85, Math.floor((fullContent.length / 2000) * 100));
                  if (estimatedProgress < 10) estimatedProgress = Math.max(10, chunkCount * 2);
                  
                  console.log(`[DEBUG] Sending progress update: ${estimatedProgress}% (chunk ${chunkCount}, content: ${fullContent.length} chars)`);
                  
                  progressCallback({
                    type: 'streaming',
                    progress: estimatedProgress,
                    reasoning: fullContent,
                    message: `Analyzing... (${Math.floor(fullContent.length / 100)} insights gathered)`,
                    sources: []
                  });
                  
                  lastProgressUpdate = now;
                }
              }
            } catch (parseError) {
              console.warn('[DEBUG] Failed to parse streaming chunk:', parseError);
            }
          }
        }
      }
    } finally {
      reader.releaseLock();
    }

    // Check if content appears complete
    const isComplete = fullContent.length > 100 && (
      fullContent.trim().endsWith('.') || 
      fullContent.trim().endsWith('!') || 
      fullContent.trim().endsWith('?') ||
      fullContent.includes('</think>')
    );

    console.log('[DEBUG] Final content check:', {
      contentLength: fullContent.length,
      isComplete: isComplete,
      lastChars: fullContent.slice(-100)
    });

    // Send completion update
    if (progressCallback) {
      console.log('[DEBUG] Sending completion update with full content');
      progressCallback({
        type: 'complete',
        progress: 100,
        reasoning: fullContent,
        message: 'Analysis complete!',
        sources: []
      });
    }

    const responseTime = Date.now() - startTime;
    console.log('[DEBUG] Streaming request completed:', {
      responseTime: `${responseTime}ms`,
      contentLength: fullContent.length,
      totalChunks: chunkCount,
      isComplete: isComplete
    });

    if (!isComplete) {
      console.warn('[DEBUG] ⚠️ RESEARCH CONTENT MAY BE TRUNCATED:', {
        contentLength: fullContent.length,
        lastChars: fullContent.slice(-100),
        tokensUsed: Math.ceil(fullContent.length / 4),
        maxTokens: data.max_tokens
      });
    }

    // Return a mock response structure that matches what the calling code expects
    return {
      id: `stream_${Date.now()}`,
      model: data.model,
      choices: [{
        message: {
          content: fullContent
        }
      }],
      usage: {
        total_tokens: Math.ceil(fullContent.length / 4)
      }
    };
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
    console.log('🔥 Optimized reasoning query called:', { query, contextLength: context?.length });
    
    const messages: any[] = [];
    const progressUpdates: ProgressUpdate[] = [];
    
    // Enhanced system prompt for reasoning
    const systemPrompt = `You are Monday, an AI assistant with advanced reasoning capabilities.

REASONING INSTRUCTIONS:
- Think through the problem step-by-step in <think> tags
- Show your complete reasoning process including analysis, connections, and conclusions
- After your reasoning, provide a clear, comprehensive final answer
- Your final answer should be well-structured and complete
- Always finish your thoughts completely - never cut off mid-sentence

RESPONSE STRUCTURE:
<think>
[Your complete step-by-step reasoning process here]
</think>

[Your final, refined answer here - this should be comprehensive and well-formatted]`;

    messages.push({
      role: 'system',
      content: systemPrompt
    });

    // Add context with proper role alternation
    if (context && context.length > 0) {
      const cleanedContext = ContextCleaner.validateAndCleanContext(context as ConversationEntry[]);
      const alternatedMessages = this.ensureAlternation(cleanedContext);
      messages.push(...alternatedMessages);
    }

    // Add the user query
    messages.push({
      role: 'user',
      content: query
    });

    console.log('[DEBUG] Reasoning query request:', {
      model: 'sonar-reasoning-pro',
      messageCount: messages.length,
      maxTokens: 2000, // Significantly increased for complete responses
      query: query.substring(0, 100)
    });

    try {
      const response = await this.makeStreamingRequest('/chat/completions', {
        model: 'sonar-reasoning-pro',
        messages: messages,
        max_tokens: 2000, // Increased from 600 to 2000 tokens
        temperature: 0.3,
        top_p: 0.9,
        stream: true,
        web_search_options: {
          search_context_size: "high" // Use high context for reasoning
        }
      }, progressCallback);

      // Parse the reasoning response to extract both thinking and final answer
      const parsedResponse = this.parseReasoningResponse(response.content);
      
      console.log('[DEBUG] Reasoning response parsed:', {
        hasThinking: !!parsedResponse.thinking,
        hasFinalAnswer: !!parsedResponse.finalAnswer,
        thinkingLength: parsedResponse.thinking?.length || 0,
        finalAnswerLength: parsedResponse.finalAnswer?.length || 0,
        totalLength: response.content.length
      });

      return {
        id: `reasoning_${Date.now()}`,
        model: 'sonar-reasoning-pro',
        content: parsedResponse.finalAnswer || response.content, // Use final answer for TTS
        fullContent: response.content, // Keep full content for display
        reasoning: parsedResponse.thinking ? this.extractReasoningSteps(parsedResponse.thinking) : undefined, // Convert to ReasoningStep[]
        sources: response.sources || [],
        metadata: {
          tokensUsed: response.usage?.total_tokens || 0,
          responseTime: 0,
          isThinking: true,
          isStreaming: true,
          ttsMessage: parsedResponse.finalAnswer ? 'Analysis complete! Here are my findings.' : 'I\'m thinking through this step by step.'
        }
      };
    } catch (error) {
      console.error('[ERROR] Reasoning query failed:', error);
      throw error;
    }
  }

  // New method to parse reasoning responses
  private parseReasoningResponse(content: string): { thinking?: string; finalAnswer?: string } {
    const thinkRegex = /<think>([\s\S]*?)<\/think>/;
    const thinkMatch = content.match(thinkRegex);
    
    let thinking: string | undefined;
    let finalAnswer: string | undefined;
    
    if (thinkMatch) {
      thinking = thinkMatch[1].trim();
      // Extract content after </think> tag as final answer
      const afterThink = content.split('</think>')[1];
      if (afterThink && afterThink.trim()) {
        finalAnswer = afterThink.trim();
      }
    }
    
    // If no think tags found, treat entire content as final answer
    if (!thinking && !finalAnswer) {
      finalAnswer = content;
    }
    
    return { thinking, finalAnswer };
  }

  public async deepResearch(query: string, context?: string[] | ConversationEntry[], progressCallback?: (update: ProgressUpdate) => void): Promise<PerplexityResponse> {
    console.log('🔥 Optimized deep research called:', { query, contextLength: context?.length });
    
    const messages: any[] = [];
    const progressUpdates: ProgressUpdate[] = [];
    
    // Enhanced system prompt for research
    const systemPrompt = `You are Monday, conducting comprehensive research analysis.

RESEARCH INSTRUCTIONS:
- Conduct thorough research using multiple high-quality sources
- Synthesize information from diverse perspectives
- Provide comprehensive analysis with proper citations
- Structure your response clearly with sections and key findings
- Always complete your analysis - never cut off mid-sentence

RESPONSE STRUCTURE:
- Opening: Brief context and research scope
- Main findings: Key insights with source backing
- Analysis: Synthesis of information and implications
- Conclusion: Summary and recommendations for further exploration

Focus on accuracy, depth, and practical insights.`;

    messages.push({
      role: 'system',
      content: systemPrompt
    });

    // Add context with proper role alternation
    if (context && context.length > 0) {
      const cleanedContext = ContextCleaner.validateAndCleanContext(context as ConversationEntry[]);
      const alternatedMessages = this.ensureAlternation(cleanedContext);
      messages.push(...alternatedMessages);
    }

    // Add the user query
    messages.push({
      role: 'user',
      content: query
    });

    console.log('[DEBUG] Deep research request:', {
      model: 'sonar-deep-research',
      messageCount: messages.length,
      maxTokens: 2500, // Increased for comprehensive research
      query: query.substring(0, 100)
    });

    try {
      const response = await this.makeStreamingRequest('/chat/completions', {
        model: 'sonar-deep-research',
        messages: messages,
        max_tokens: 2500, // Increased from 700 to 2500 tokens for comprehensive research
        temperature: 0.2,
        top_p: 0.9,
        stream: true,
        web_search_options: {
          search_context_size: "high" // Use high context for deep research
        }
      }, progressCallback);

      console.log('[DEBUG] Deep research response received:', {
        contentLength: response.content?.length || 0,
        hasContent: !!response.content,
        totalTokens: response.usage?.total_tokens || 0
      });

      return {
        id: `research_${Date.now()}`,
        model: 'sonar-deep-research',
        content: response.content || 'Research completed. Please check the full analysis.',
        fullContent: response.content,
        sources: response.sources || [],
        metadata: {
          tokensUsed: response.usage?.total_tokens || 0,
          responseTime: 0,
          isResearching: true,
          isStreaming: true,
          ttsMessage: 'I\'ve completed comprehensive research on this topic. Here\'s what I found.'
        }
      };
    } catch (error) {
      console.error('[ERROR] Deep research failed:', error);
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