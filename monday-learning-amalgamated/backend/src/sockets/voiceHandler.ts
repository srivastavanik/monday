import { perplexityService } from '../services/perplexity.js'
import { logger } from '../utils/logger.js'

// Monday command patterns
const MONDAY_PATTERNS = {
  greeting: /^monday,?\s+(hello|hi|hey|start)/i,
  basic: /^monday,?\s+(.+)/i,
  reasoning: /^monday,?\s+(think\s+about|reason\s+through|analyze)\s+(.+)/i,
  deepResearch: /^monday,?\s+(deep\s+dive|research|investigate)\s+(.+)/i,
  spatial: /^monday,?\s+(bring\s+this\s+closer|push\s+that\s+away|move\s+.+)/i,
  focus: /^monday,?\s+focus\s+mode/i,
  context: /^monday,?\s+(what\s+about|tell\s+me\s+about)\s+(.+)\s+from\s+earlier/i
}

function parseCommand(text: string) {
  const cleanText = text.trim().toLowerCase()
  
  // Check for Monday activation
  if (!cleanText.startsWith('monday')) {
    return null
  }

  let commandType = 'basic'
  let content = ''

  if (MONDAY_PATTERNS.greeting.test(cleanText)) {
    commandType = 'greeting'
    content = 'greeting'
  } else if (MONDAY_PATTERNS.reasoning.test(cleanText)) {
    commandType = 'reasoning'
    const match = cleanText.match(MONDAY_PATTERNS.reasoning)
    content = match?.[2] || ''
  } else if (MONDAY_PATTERNS.deepResearch.test(cleanText)) {
    commandType = 'deepResearch'
    const match = cleanText.match(MONDAY_PATTERNS.deepResearch)
    content = match?.[2] || ''
  } else if (MONDAY_PATTERNS.spatial.test(cleanText)) {
    commandType = 'spatial'
    content = cleanText.replace(/^monday,?\s+/, '')
  } else if (MONDAY_PATTERNS.focus.test(cleanText)) {
    commandType = 'focus'
    content = 'focus'
  } else if (MONDAY_PATTERNS.context.test(cleanText)) {
    commandType = 'context'
    const match = cleanText.match(MONDAY_PATTERNS.context)
    content = match?.[2] || ''
  } else if (MONDAY_PATTERNS.basic.test(cleanText)) {
    commandType = 'basic'
    const match = cleanText.match(MONDAY_PATTERNS.basic)
    content = match?.[1] || ''
  }

  return {
    type: commandType,
    content: content.trim(),
    originalText: text,
    timestamp: Date.now()
  }
}

export const handleVoiceQuery = (socket: any, io: any) => {
  socket.on('voice_command', async (data: any) => {
    try {
      logger.info('Voice command received', { 
        socketId: socket.id, 
        command: data.command?.substring(0, 50) 
      })
      
      const command = parseCommand(data.command || '')
      
      if (!command) {
        socket.emit('monday_response', {
          type: 'error',
          content: 'Please start your command with "Monday"',
          timestamp: Date.now()
        })
        return
      }

      // Handle different command types
      switch (command.type) {
        case 'greeting':
          socket.emit('monday_response', {
            type: 'greeting',
            content: "Hello! I'm Monday, your AI learning companion. I'm here to help you explore knowledge in three dimensions. What would you like to learn about today?",
            timestamp: Date.now()
          })
          break

        case 'basic':
          if (command.content) {
            const response = await perplexityService.processQuery({
              query: command.content,
              mode: 'basic',
              sessionId: data.sessionId
            })
            
            socket.emit('monday_response', {
              type: 'basic_response',
              content: response.content,
              citations: response.citations,
              metadata: response.metadata,
              timestamp: Date.now()
            })
          }
          break

        case 'reasoning':
          if (command.content) {
            const response = await perplexityService.processQuery({
              query: command.content,
              mode: 'reasoning',
              sessionId: data.sessionId
            })
            
            socket.emit('monday_response', {
              type: 'reasoning_response',
              content: response.content,
              reasoning: response.reasoning,
              citations: response.citations,
              metadata: response.metadata,
              timestamp: Date.now()
            })
          }
          break

        case 'deepResearch':
          if (command.content) {
            const response = await perplexityService.processQuery({
              query: command.content,
              mode: 'research',
              sessionId: data.sessionId
            })
            
            socket.emit('monday_response', {
              type: 'research_response',
              content: response.content,
              sources: response.sources,
              citations: response.citations,
              metadata: response.metadata,
              timestamp: Date.now()
            })
          }
          break

        case 'spatial':
          socket.emit('spatial_command', {
            type: 'spatial_manipulation',
            command: command.content,
            timestamp: Date.now()
          })
          break

        case 'focus':
          socket.emit('focus_mode', {
            enabled: true,
            timestamp: Date.now()
          })
          break

        case 'context':
          socket.emit('monday_response', {
            type: 'context_response',
            content: `I understand you're asking about "${command.content}" from our earlier discussion. Let me recall that for you.`,
            timestamp: Date.now()
          })
          break

        default:
          socket.emit('monday_response', {
            type: 'general_response',
            content: `I heard you say: "${command.originalText}". How can I help you learn today?`,
            timestamp: Date.now()
          })
      }

    } catch (error) {
      logger.error('Voice command processing failed:', error)
      socket.emit('monday_response', {
        type: 'error',
        content: 'I encountered an error processing your request. Please try again.',
        timestamp: Date.now()
      })
    }
  })

  // Handle typing indicator
  socket.on('monday_thinking', (data: any) => {
    socket.emit('monday_status', {
      type: 'processing',
      message: 'Monday is thinking...',
      timestamp: Date.now()
    })
  })
} 