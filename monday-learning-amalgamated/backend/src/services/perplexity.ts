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

class PerplexityService {
  private apiKey: string
  private baseUrl = 'https://api.perplexity.ai'
  
  constructor() {
    this.apiKey = process.env.PERPLEXITY_API_KEY || ''
    if (!this.apiKey) {
      throw new Error('PERPLEXITY_API_KEY environment variable is required')
    }
  }

  private async makeRequest(endpoint: string, data: any): Promise<any> {
    const startTime = Date.now()
    
    try {
      const response = await axios.post(`${this.baseUrl}${endpoint}`, data, {
        headers: {
          'Authorization': `Bearer ${this.apiKey}`,
          'Content-Type': 'application/json'
        },
        timeout: 30000 // 30 second timeout
      })
      
      const responseTime = Date.now() - startTime
      logger.info('Perplexity API request completed', {
        endpoint,
        responseTime: `${responseTime}ms`,
        statusCode: response.status
      })
      
      return response.data
    } catch (error: any) {
      const responseTime = Date.now() - startTime
      logger.error('Perplexity API request failed', {
        endpoint,
        responseTime: `${responseTime}ms`,
        error: error.message,
        status: error.response?.status,
        responseData: error.response?.data ? JSON.stringify(error.response.data) : 'No response data'
      })
      throw error
    }
  }

  async basicQuery(query: string, context?: string[]): Promise<PerplexityResponse> {
    const messages = [
      {
        role: 'system',
        content: `You are Monday, an advanced AI learning companion for VR education, powered by Perplexity Sonar.

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
- Always end with engagement (questions, suggestions, or offers to explore more)`
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

    // Add the user query
    messages.push({
      role: 'user',
      content: query
    })

    const requestData = {
      model: 'sonar-pro',
      messages: messages,
      max_tokens: 300,
      temperature: 0.3
      // Removed search_domain_filter as it was causing errors with invalid domain names
      // The API works best with default search settings according to documentation
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    return {
      id: result.id || 'basic_query',
      model: result.model || 'sonar-pro',
      content: result.choices?.[0]?.message?.content || 'No response generated',
      citations: this.extractCitations(result),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0 // Will be filled by caller
      }
    }
  }

  async reasoningQuery(query: string, context?: string[]): Promise<PerplexityResponse> {
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
      temperature: 0.2
      // Using default search settings for best results
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    return {
      id: result.id || 'reasoning_query',
      model: result.model || 'sonar-reasoning-pro',
      content: result.choices?.[0]?.message?.content || 'No response generated',
      citations: this.extractCitations(result),
      reasoning: this.extractReasoningSteps(result.choices?.[0]?.message?.content || ''),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0
      }
    }
  }

  async deepResearch(query: string): Promise<PerplexityResponse> {
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
      temperature: 0.3
      // Deep research model works best with minimal configuration
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    return {
      id: result.id || 'research_query',
      model: result.model || 'sonar-deep-research',
      content: result.choices?.[0]?.message?.content || 'No response generated',
      citations: this.extractCitations(result),
      sources: this.extractSources(result),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0
      }
    }
  }

  private extractCitations(result: any): Citation[] {
    const citations: Citation[] = []
    
    // Perplexity API now automatically includes citations in responses
    // Check for citations in the response data
    if (result.citations) {
      result.citations.forEach((citation: any, index: number) => {
        citations.push({
          id: `citation_${index}`,
          url: citation.url || '',
          title: citation.title || 'Untitled',
          snippet: citation.text || citation.snippet || '',
          publishedDate: citation.published_date,
          domain: this.extractDomain(citation.url || '')
        })
      })
    }
    
    return citations
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

  async processQuery(queryData: PerplexityQuery): Promise<PerplexityResponse> {
    const startTime = Date.now()
    
    let response: PerplexityResponse
    
    switch (queryData.mode) {
      case 'basic':
        response = await this.basicQuery(queryData.query, queryData.context)
        break
      case 'reasoning':
        response = await this.reasoningQuery(queryData.query, queryData.context)
        break
      case 'research':
        response = await this.deepResearch(queryData.query)
        break
      default:
        throw new Error(`Unsupported query mode: ${queryData.mode}`)
    }
    
    // Add response time to metadata
    response.metadata.responseTime = Date.now() - startTime
    
    // Log the query for analytics
    logger.info('Perplexity query processed', {
      mode: queryData.mode,
      query: queryData.query.substring(0, 100) + '...',
      responseTime: response.metadata.responseTime,
      tokensUsed: response.metadata.tokensUsed,
      citationCount: response.citations?.length || 0,
      sessionId: queryData.sessionId
    })
    
    return response
  }
}

export const perplexityService = new PerplexityService() 