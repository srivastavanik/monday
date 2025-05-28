export interface APIResponse {
  id: string
  mode: "basic" | "reasoning" | "deep-research"
  query: string
  content: string
  reasoning_steps?: ReasoningStep[]
  sources?: Source[]
  confidence?: number
  timestamp: Date
}

export interface ReasoningStep {
  id: string
  step: string
  reasoning: string
  confidence: number
}

export interface Source {
  id: string
  title: string
  url: string
  relevance: number
  excerpt: string
}

export interface LearningNode {
  id: string
  topic: string
  understanding: number
  connections: string[]
  position: { x: number; y: number }
}
