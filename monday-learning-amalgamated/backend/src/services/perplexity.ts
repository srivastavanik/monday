import axios from 'axios'
import { logger } from '../utils/logger.js'

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

  public async reasoningQuery(query: string, context?: string[]): Promise<PerplexityResponse> {
    const messages = [
      {
        role: 'system',
        content: `You are Monday, an AI learning companion with advanced reasoning capabilities.

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
- Make complex topics approachable`
      }
    ]

    // Add context messages if provided
    if (context && context.length > 0) {
      context.forEach(ctx => {
        messages.push({
          role: 'assistant',
          content: ctx
        })
      })
    }

    // Add the user query with reasoning prompt
    messages.push({
      role: 'user',
      content: `Please think through this step by step: ${query}`
    })

    const requestData = {
      model: 'sonar-reasoning-pro',
      messages: messages,
      max_tokens: 500,
      temperature: 0.2,
      stream: false
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    const fullContent = result.choices?.[0]?.message?.content || 'No response generated'
    const sentences = fullContent.split(/[.!?]+/).filter(s => s.trim().length > 0)
    const ttsContent = sentences.slice(0, 2).join('. ').trim()
    const finalTtsContent = ttsContent.endsWith('.') ? ttsContent : ttsContent + '.'
    
    return {
      id: result.id || 'reasoning_query',
      model: result.model || 'sonar-reasoning-pro',
      content: finalTtsContent,
      fullContent: fullContent,
      citations: this.extractCitations(result),
      reasoning: this.extractReasoningSteps(fullContent),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0
      }
    }
  }

  public async deepResearch(query: string): Promise<PerplexityResponse> {
    const requestData = {
      model: 'sonar-deep-research',
      messages: [
        {
          role: 'system',
          content: `You are Monday, conducting comprehensive research analysis.

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
- Maintain conversational tone despite depth`
        },
        {
          role: 'user',
          content: `Conduct a comprehensive research analysis on: ${query}`
        }
      ],
      max_tokens: 800,
      temperature: 0.3,
      stream: false
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    const fullContent = result.choices?.[0]?.message?.content || 'No response generated'
    const sentences = fullContent.split(/[.!?]+/).filter(s => s.trim().length > 0)
    const ttsContent = sentences.slice(0, 2).join('. ').trim()
    const finalTtsContent = ttsContent.endsWith('.') ? ttsContent : ttsContent + '.'
    
    return {
      id: result.id || 'research_query',
      model: result.model || 'sonar-deep-research',
      content: finalTtsContent,
      fullContent: fullContent,
      citations: this.extractCitations(result),
      sources: this.extractSources(result),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0
      }
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
    // Extract first 1-2 sentences for TTS
    const sentences = fullContent.split(/[.!?]+/).filter(s => s.trim().length > 0)
    let shortResponse = sentences.slice(0, 2).join('. ').trim()
    
    // Ensure it ends with punctuation
    if (!shortResponse.endsWith('.') && !shortResponse.endsWith('!') && !shortResponse.endsWith('?')) {
      shortResponse += '.'
    }
    
    // If too long, truncate to ~150 characters
    if (shortResponse.length > 150) {
      shortResponse = shortResponse.substring(0, 147) + '...'
    }
    
    return shortResponse
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