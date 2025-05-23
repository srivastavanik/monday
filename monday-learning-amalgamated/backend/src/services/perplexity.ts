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
        status: error.response?.status
      })
      throw error
    }
  }

  async basicQuery(query: string, context?: string[]): Promise<PerplexityResponse> {
    const requestData = {
      model: 'llama-3.1-sonar-small-128k-online',
      messages: [
        {
          role: 'system',
          content: 'You are Monday, a helpful AI learning companion. Provide clear, educational responses with proper citations.'
        },
        ...context ? context.map(ctx => ({
          role: 'assistant',
          content: ctx
        })) : [],
        {
          role: 'user',
          content: query
        }
      ],
      return_citations: true,
      search_domain_filter: ['educational', 'academic'],
      temperature: 0.2
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    return {
      id: result.id,
      model: result.model,
      content: result.choices[0].message.content,
      citations: this.extractCitations(result),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0 // Will be filled by caller
      }
    }
  }

  async reasoningQuery(query: string, context?: string[]): Promise<PerplexityResponse> {
    const requestData = {
      model: 'llama-3.1-sonar-large-128k-online',
      messages: [
        {
          role: 'system',
          content: `You are Monday, an AI learning companion with advanced reasoning capabilities. 
                   Break down complex problems step by step, show your reasoning process, and provide confidence scores for each step.
                   Format your response with clear reasoning steps.`
        },
        ...context ? context.map(ctx => ({
          role: 'assistant',
          content: ctx
        })) : [],
        {
          role: 'user',
          content: `Please think through this step by step: ${query}`
        }
      ],
      return_citations: true,
      search_domain_filter: ['educational', 'academic'],
      temperature: 0.1,
      return_related_questions: true
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    return {
      id: result.id,
      model: result.model,
      content: result.choices[0].message.content,
      citations: this.extractCitations(result),
      reasoning: this.extractReasoningSteps(result.choices[0].message.content),
      metadata: {
        tokensUsed: result.usage?.total_tokens || 0,
        responseTime: 0
      }
    }
  }

  async deepResearch(query: string): Promise<PerplexityResponse> {
    const requestData = {
      model: 'llama-3.1-sonar-huge-128k-online',
      messages: [
        {
          role: 'system',
          content: `You are Monday, conducting deep research analysis. Provide comprehensive, multi-source analysis with:
                   1. Multiple perspectives on the topic
                   2. Detailed source analysis
                   3. Connections between different sources
                   4. Critical evaluation of information
                   5. Synthesis of findings`
        },
        {
          role: 'user',
          content: `Conduct a deep research analysis on: ${query}`
        }
      ],
      return_citations: true,
      search_recency_filter: 'month',
      search_domain_filter: ['educational', 'academic', 'news'],
      temperature: 0.3,
      max_tokens: 4000
    }

    const result = await this.makeRequest('/chat/completions', requestData)
    
    return {
      id: result.id,
      model: result.model,
      content: result.choices[0].message.content,
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
    
    if (result.citations) {
      result.citations.forEach((citation: any, index: number) => {
        citations.push({
          id: `citation_${index}`,
          url: citation.url || '',
          title: citation.title || 'Untitled',
          snippet: citation.text || '',
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
          snippet: citation.text || '',
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