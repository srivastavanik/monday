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
}

export interface PerplexityResponse {
  id: string
  model: string
  content: string
  fullContent?: string
  citations?: Citation[]
  reasoning?: ReasoningStep[]
  sources?: Source[]
  metadata: {
    tokensUsed: number
    responseTime: number
    confidence?: number
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
        timeout: 30000
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
        max_tokens: 300,
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

  public async reasoningQuery(query: string, context?: string[] | ConversationEntry[]): Promise<PerplexityResponse> {
    console.log('🔥 NEW REASONING QUERY METHOD CALLED!', { query, contextLength: context?.length });
    
    const messages: any[] = [];
    
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
- Keep total response under 500 tokens for voice delivery

Educational Focus:
- Help users understand not just what, but why and how
- Encourage critical thinking
- Make complex topics approachable`;

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
        // Convert old string format to new format
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
      max_tokens: 500,
      temperature: 0.2,
      stream: false
    };

    try {
      const result = await this.makeRequest('/chat/completions', requestData);
      
      const fullContent = result.choices?.[0]?.message?.content || 'No response generated';
      const shortResponse = this.createShortTTSResponse(fullContent, query);
      
      return {
        id: result.id || 'reasoning_query',
        model: result.model || 'sonar-reasoning-pro',
        content: shortResponse,
        fullContent: fullContent,
        citations: this.extractCitations(result),
        reasoning: this.extractReasoningSteps(fullContent),
        metadata: {
          tokensUsed: result.usage?.total_tokens || 0,
          responseTime: 0
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

  public async deepResearch(query: string, context?: string[] | ConversationEntry[]): Promise<PerplexityResponse> {
    const messages: any[] = [];
    
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

Voice-Friendly Delivery:
- Use clear, flowing language suitable for TTS
- Break up long sections with natural pauses
- Avoid excessive technical jargon without explanation
- Maintain conversational tone despite depth`;

    messages.push({
      role: 'system',
      content: systemPrompt
    });

    // Process context with new structure
    if (context && Array.isArray(context)) {
      let structuredContext: ConversationEntry[];
      
      // Handle both old string format and new ConversationEntry format
      if (context.length > 0 && typeof context[0] === 'string') {
        // Convert old string format to new format
        structuredContext = ContextCleaner.convertLegacyContext(context as string[]);
      } else {
        structuredContext = context as ConversationEntry[];
      }
      
      // Clean the context
      const cleanedMessages = ContextCleaner.validateAndCleanContext(structuredContext);
      
      // Ensure perfect alternation
      const alternatingMessages = this.ensureAlternation(cleanedMessages);
      
      // Add to messages
      messages.push(...alternatingMessages);
    }

    // Add the user query
    messages.push({
      role: 'user',
      content: `Conduct a comprehensive research analysis on: ${query}`
    });

    // Final validation
    this.validateMessageStructure(messages);
    
    // Log for debugging
    this.logApiRequest('/chat/completions', {
      model: 'sonar-deep-research',
      messages: messages
    });

    const requestData = {
      model: 'sonar-deep-research',
      messages: messages,
      max_tokens: 800,
      temperature: 0.3,
      stream: false
    };

    try {
      const result = await this.makeRequest('/chat/completions', requestData);
      
      const fullContent = result.choices?.[0]?.message?.content || 'No response generated';
      const shortResponse = this.createShortTTSResponse(fullContent, query);
      
      return {
        id: result.id || 'research_query',
        model: result.model || 'sonar-deep-research',
        content: shortResponse,
        fullContent: fullContent,
        citations: this.extractCitations(result),
        sources: this.extractSources(result),
        metadata: {
          tokensUsed: result.usage?.total_tokens || 0,
          responseTime: 0
        }
      };
    } catch (error: any) {
      console.error('Perplexity API Error:', error.response?.data || error.message);
      throw error;
    }
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
          sources: []
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