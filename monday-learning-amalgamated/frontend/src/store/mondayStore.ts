import { create } from 'zustand'
import { devtools } from 'zustand/middleware'

// Types for the Monday store
export interface MondayPanel {
  id: string
  type: 'content' | 'reasoning' | 'research' | 'visualization' | 'citations' | 'progressive_reasoning' | 'progressive_research' | 'thinking'
  position: [number, number, number]
  rotation: [number, number, number]
  content: string
  title: string
  isActive: boolean
  opacity: number
  citations?: Citation[]
  reasoning?: any[]
  createdAt: number
  updatedAt?: number
  // Additional properties for progressive panels
  progressive?: boolean
  fullContent?: string
  sources?: any[]
  isThinking?: boolean
  isResearching?: boolean
  model?: string
}

export interface Citation {
  id: string
  url: string
  title: string
  snippet: string
  position: [number, number, number]
}

export interface ReasoningStep {
  id: string
  content: string
  confidence: number
  connections: string[]
  position: [number, number, number]
  timestamp: number
}

export interface ResearchNode {
  id: string
  title: string
  type: 'source' | 'concept' | 'connection'
  position: [number, number, number]
  connections: string[]
  metadata: Record<string, any>
}

export interface SessionState {
  isActive: boolean
  hasCompletedIntro: boolean
  currentMode: 'basic' | 'reasoning' | 'research' | 'focus'
  learningMomentum: number
  totalQueries: number
  sessionStartTime: number
  lastInteractionTime: number
  topicsExplored: string[]
}

export interface PerformanceMetrics {
  fps: number
  frameTime: number
  memoryUsage: number
  timestamp: number
}

export interface MondayState {
  // Connection state
  isConnected: boolean
  socket: any
  
  // Session management
  sessionState: SessionState
  
  // Spatial data
  panels: MondayPanel[]
  activePanel: string | null
  reasoningChain: ReasoningStep[]
  researchWeb: ResearchNode[]
  
  // Voice and interaction
  isListening: boolean
  currentTranscript: string
  lastCommand: string | null
  mondayResponse: string | null
  conversationActive: boolean
  
  // Performance monitoring
  performanceMetrics: PerformanceMetrics
  
  // Visual state
  focusMode: boolean
  spatialLayout: 'default' | 'focus' | 'research'
  
  // Learning tracking
  knowledgeConstellation: Array<{
    topic: string
    understanding: number
    position: [number, number, number]
    connections: string[]
  }>
  
  // Actions
  setConnected: (connected: boolean) => void
  setSocket: (socket: any) => void
  initializeSession: (config: any) => void
  addPanel: (panel: Partial<MondayPanel> & { type: string; content: string; title: string; position: [number, number, number]; rotation: [number, number, number] }) => void
  updatePanel: (id: string, updates: Partial<MondayPanel>) => void
  removePanel: (id: string) => void
  setActivePanel: (id: string | null) => void
  addReasoningStep: (step: Omit<ReasoningStep, 'id' | 'timestamp'>) => void
  clearReasoningChain: () => void
  addResearchNode: (node: Omit<ResearchNode, 'id'>) => void
  clearResearchWeb: () => void
  setTranscript: (transcript: string) => void
  setMondayResponse: (response: string) => void
  setLastCommand: (command: string) => void
  setListening: (listening: boolean) => void
  setPerformanceMetrics: (metrics: PerformanceMetrics) => void
  setFocusMode: (enabled: boolean) => void
  setSpatialLayout: (layout: 'default' | 'focus' | 'research') => void
  updateLearningMomentum: (change: number) => void
  addToKnowledgeConstellation: (topic: string, understanding: number) => void
  setConversationActive: (active: boolean) => void
}

export const useMondayStore = create<MondayState>()(
  devtools(
    (set, get) => ({
      // Initial state
      isConnected: false,
      socket: null,
      
      sessionState: {
        isActive: false,
        hasCompletedIntro: false,
        currentMode: 'basic',
        learningMomentum: 0,
        totalQueries: 0,
        sessionStartTime: Date.now(),
        lastInteractionTime: Date.now(),
        topicsExplored: []
      },
      
      panels: [],
      activePanel: null,
      reasoningChain: [],
      researchWeb: [],
      
      isListening: false,
      currentTranscript: '',
      lastCommand: null,
      mondayResponse: null,
      conversationActive: false,
      
      performanceMetrics: {
        fps: 60,
        frameTime: 16.67,
        memoryUsage: 0,
        timestamp: Date.now()
      },
      
      focusMode: false,
      spatialLayout: 'default',
      
      knowledgeConstellation: [],
      
      // Actions
      setConnected: (connected) => set({ isConnected: connected }),
      
      setSocket: (socket) => set({ socket }),
      
      setConversationActive: (active) => set({ conversationActive: active }),
      
      initializeSession: (config) => set((state) => ({
        sessionState: {
          ...state.sessionState,
          isActive: true,
          sessionStartTime: Date.now(),
          lastInteractionTime: Date.now()
        }
      })),
      
      addPanel: (panelData) => {
        console.log('Store: Adding panel with data:', panelData)
        
        const newPanel: MondayPanel = {
          // Set defaults first
          id: `panel_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
          createdAt: Date.now(),
          updatedAt: Date.now(),
          isActive: false,
          opacity: 1,
          // Then spread the incoming data (which may override defaults)
          ...panelData,
          // Ensure required fields are present
          type: panelData.type as any, // Cast to handle 'citations' type
        }
        
        console.log('Store: Created panel object:', newPanel)
        
        set((state) => {
          const newPanels = [...state.panels, newPanel]
          console.log('Store: Updated panels array length:', newPanels.length)
          return {
            panels: newPanels,
            activePanel: newPanel.isActive ? newPanel.id : state.activePanel
          }
        })
      },
      
      updatePanel: (id, updates) => set((state) => ({
        panels: state.panels.map(panel => 
          panel.id === id 
            ? { ...panel, ...updates, updatedAt: Date.now() }
            : panel
        )
      })),
      
      removePanel: (id) => set((state) => ({
        panels: state.panels.filter(panel => panel.id !== id),
        activePanel: state.activePanel === id ? null : state.activePanel
      })),
      
      setActivePanel: (id) => set({ activePanel: id }),
      
      addReasoningStep: (stepData) => {
        const newStep: ReasoningStep = {
          ...stepData,
          id: `step_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
          timestamp: Date.now()
        }
        
        set((state) => ({
          reasoningChain: [...state.reasoningChain, newStep]
        }))
      },
      
      clearReasoningChain: () => set({ reasoningChain: [] }),
      
      addResearchNode: (nodeData) => {
        const newNode: ResearchNode = {
          ...nodeData,
          id: `node_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`
        }
        
        set((state) => ({
          researchWeb: [...state.researchWeb, newNode]
        }))
      },
      
      clearResearchWeb: () => set({ researchWeb: [] }),
      
      setTranscript: (transcript) => set({ currentTranscript: transcript }),
      
      setMondayResponse: (response) => set({ mondayResponse: response }),
      
      setLastCommand: (command) => set({ lastCommand: command }),
      
      setListening: (listening) => set({ isListening: listening }),
      
      setPerformanceMetrics: (metrics) => set({ performanceMetrics: metrics }),
      
      setFocusMode: (enabled) => {
        set((state) => ({
          focusMode: enabled,
          spatialLayout: enabled ? 'focus' : 'default',
          panels: state.panels.map(panel => ({
            ...panel,
            opacity: enabled ? (panel.isActive ? 1 : 0.15) : 1
          }))
        }))
      },
      
      setSpatialLayout: (layout) => set({ spatialLayout: layout }),
      
      updateLearningMomentum: (change) => set((state) => ({
        sessionState: {
          ...state.sessionState,
          learningMomentum: Math.max(0, Math.min(100, state.sessionState.learningMomentum + change)),
          lastInteractionTime: Date.now()
        }
      })),
      
      addToKnowledgeConstellation: (topic, understanding) => {
        const existingIndex = get().knowledgeConstellation.findIndex(
          item => item.topic.toLowerCase() === topic.toLowerCase()
        )
        
        if (existingIndex >= 0) {
          // Update existing topic
          set((state) => ({
            knowledgeConstellation: state.knowledgeConstellation.map((item, index) =>
              index === existingIndex
                ? { ...item, understanding: Math.max(item.understanding, understanding) }
                : item
            )
          }))
        } else {
          // Add new topic to constellation
          const position: [number, number, number] = [
            (Math.random() - 0.5) * 10,
            Math.random() * 2 + 1,
            (Math.random() - 0.5) * 10
          ]
          
          set((state) => ({
            knowledgeConstellation: [
              ...state.knowledgeConstellation,
              {
                topic,
                understanding,
                position,
                connections: []
              }
            ]
          }))
        }
      }
    }),
    {
      name: 'monday-store',
      version: 1
    }
  )
) 