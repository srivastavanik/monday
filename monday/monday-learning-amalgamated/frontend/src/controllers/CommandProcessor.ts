interface CommandEvent {
  id: string
  transcript: string
  timestamp: number
  processed: boolean
}

interface ConversationContext {
  active: boolean
  startTime: number | null
  lastCommandTime: number | null
  commandCount: number
  conversationId: string | null
}

type ConversationListener = (context: ConversationContext) => void

class CommandProcessor {
  private static instance: CommandProcessor
  private commandQueue: CommandEvent[] = []
  private processing = false
  private conversationContext: ConversationContext = {
    active: false,
    startTime: null,
    lastCommandTime: null,
    commandCount: 0,
    conversationId: null
  }
  
  // Event emitter for UI updates only
  private listeners: Set<ConversationListener> = new Set()
  
  public static getInstance(): CommandProcessor {
    if (!CommandProcessor.instance) {
      CommandProcessor.instance = new CommandProcessor()
    }
    return CommandProcessor.instance
  }
  
  // Called by voice recognition - no React state checks
  public queueCommand(transcript: string, timestamp: number): void {
    const event: CommandEvent = {
      id: this.generateId(),
      transcript,
      timestamp,
      processed: false
    }
    
    console.log(`🎯 CommandProcessor: Queued command: "${transcript}"`)
    this.commandQueue.push(event)
    this.processQueue()
  }
  
  private async processQueue(): Promise<void> {
    if (this.processing || this.commandQueue.length === 0) return
    
    this.processing = true
    
    while (this.commandQueue.length > 0) {
      const event = this.commandQueue.shift()!
      await this.processCommand(event)
    }
    
    this.processing = false
  }
  
  private async processCommand(event: CommandEvent): Promise<void> {
    const normalizedTranscript = event.transcript.toLowerCase().trim()
    
    // Decision logic happens HERE, not in React
    const isActivation = normalizedTranscript.includes('hey monday')
    const isWithinConversation = this.isConversationActive()
    
    console.log(`🔍 CommandProcessor: Evaluating command: "${event.transcript}"`, {
      isActivation,
      isWithinConversation,
      conversationActive: this.conversationContext.active,
      timeSinceLastCommand: this.conversationContext.lastCommandTime ? 
        Date.now() - this.conversationContext.lastCommandTime : 'N/A'
    })
    
    if (isActivation || isWithinConversation) {
      console.log(`✅ CommandProcessor: Processing command: "${event.transcript}"`)
      
      // Update context
      if (isActivation && !this.conversationContext.active) {
        this.startConversation()
      }
      
      this.conversationContext.lastCommandTime = event.timestamp
      this.conversationContext.commandCount++
      
      // Send to backend
      await this.sendToBackend(event.transcript, isActivation)
      
      // Notify UI listeners
      this.notifyListeners()
    } else {
      console.log(`🚫 CommandProcessor: Ignoring non-conversation command: "${event.transcript}"`)
    }
    
    event.processed = true
  }
  
  public isConversationActive(): boolean {
    if (!this.conversationContext.active) return false
    
    // Conversation timeout: 5 minutes of inactivity
    const timeSinceLastCommand = Date.now() - (this.conversationContext.lastCommandTime || 0)
    if (timeSinceLastCommand > 300000) {
      console.log(`⏰ CommandProcessor: Conversation timed out after ${timeSinceLastCommand}ms`)
      this.endConversation()
      return false
    }
    
    return true
  }
  
  private startConversation(): void {
    const conversationId = this.generateId()
    console.log(`💬 CommandProcessor: Starting conversation ${conversationId}`)
    
    this.conversationContext = {
      active: true,
      startTime: Date.now(),
      lastCommandTime: Date.now(),
      commandCount: 1,
      conversationId
    }
  }
  
  public endConversation(): void {
    console.log(`🔚 CommandProcessor: Ending conversation ${this.conversationContext.conversationId}`)
    this.conversationContext.active = false
    this.notifyListeners()
  }
  
  public setConversationActive(active: boolean): void {
    if (active && !this.conversationContext.active) {
      this.startConversation()
    } else if (!active && this.conversationContext.active) {
      this.endConversation()
    }
  }
  
  private async sendToBackend(transcript: string, isExplicitTrigger: boolean): Promise<void> {
    return new Promise((resolve) => {
      try {
        const socket = (window as any).socket
        
        if (!socket || !socket.connected) {
          console.error('❌ CommandProcessor: No socket connection available')
          resolve()
          return
        }
        
        console.log(`📤 CommandProcessor: Sending to backend: "${transcript}"`)
        
        socket.emit('voice_command', {
          command: transcript,
          timestamp: Date.now(),
          isExplicitTrigger,
          conversationActive: this.conversationContext.active,
          conversationId: this.conversationContext.conversationId
        })
        
        console.log(`✅ CommandProcessor: Command sent to backend successfully`)
        
        // Don't wait for response to resolve
        setTimeout(resolve, 100)
      } catch (error) {
        console.error('❌ CommandProcessor: Error sending to backend:', error)
        resolve()
      }
    })
  }
  
  // For UI updates only
  public subscribe(listener: ConversationListener): () => void {
    this.listeners.add(listener)
    listener({ ...this.conversationContext }) // Initial state
    
    return () => {
      this.listeners.delete(listener)
    }
  }
  
  private notifyListeners(): void {
    this.listeners.forEach(listener => listener({ ...this.conversationContext }))
  }
  
  public getConversationContext(): ConversationContext {
    return { ...this.conversationContext }
  }
  
  private generateId(): string {
    return Math.random().toString(36).substring(2, 15) + Math.random().toString(36).substring(2, 15)
  }
  
  // Debug methods
  public getQueueStatus(): { queueLength: number, processing: boolean } {
    return {
      queueLength: this.commandQueue.length,
      processing: this.processing
    }
  }
  
  public clearQueue(): void {
    this.commandQueue = []
    this.processing = false
    console.log('🧹 CommandProcessor: Queue cleared')
  }
}

export { CommandProcessor }
export type { CommandEvent, ConversationContext, ConversationListener } 