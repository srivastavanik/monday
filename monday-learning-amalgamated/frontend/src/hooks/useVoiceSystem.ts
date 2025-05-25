import { useState, useEffect, useCallback, useRef } from 'react'
import { VoiceSystemController, SystemState, SystemStatus } from '../controllers/VoiceSystemController'

interface UseVoiceSystemOptions {
  onCommand?: (command: string) => void
  onError?: (error: string) => void
}

interface UseVoiceSystemReturn {
  // State
  systemState: SystemState
  isListening: boolean
  isSpeaking: boolean
  isPlaying: boolean
  transcript: string
  conversationActive: boolean
  error: string | null
  systemStatus: SystemStatus
  
  // Actions
  startConversation: () => Promise<void>
  handleCommand: (command: string) => Promise<void>
  handleTTSResponse: (text: string) => Promise<void>
  emergencyReset: () => Promise<void>
  
  // Manual overrides (for fallback UI)
  forceStartListening: () => Promise<void>
  forceStop: () => void
}

export const useVoiceSystem = (options: UseVoiceSystemOptions = {}): UseVoiceSystemReturn => {
  const controllerRef = useRef<VoiceSystemController | null>(null)
  
  // React state that mirrors controller state
  const [systemState, setSystemState] = useState<SystemState>(SystemState.IDLE)
  const [transcript, setTranscript] = useState('')
  const [conversationActive, setConversationActive] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [systemStatus, setSystemStatus] = useState<SystemStatus>({
    recognitionActive: false,
    ttsPlaying: false,
    ttsGenerating: false,
    lastTransition: Date.now(),
    errorCount: 0,
    lastError: null
  })

  // Initialize controller
  useEffect(() => {
    console.log('VoiceSystem: 🏗️ Initializing controller...')
    
    const controller = new VoiceSystemController({
      stability: 0.6,
      similarityBoost: 0.8,
      speed: 1.1
    })

    // Set up callbacks
    controller.setCallbacks({
      onStateChange: (state: SystemState) => {
        console.log('VoiceSystem: 📊 State change:', state)
        setSystemState(state)
      },
      onTranscriptChange: (newTranscript: string) => {
        console.log('VoiceSystem: 🎤 Transcript:', newTranscript)
        setTranscript(newTranscript)
        // Command processing is now handled in App.tsx via useEffect to avoid stale closures
      },
      onError: (errorMsg: string) => {
        console.error('VoiceSystem: ❌ Error:', errorMsg)
        setError(errorMsg)
        if (options.onError) {
          options.onError(errorMsg)
        }
      },
      onConversationChange: (active: boolean) => {
        console.log('VoiceSystem: 💬 Conversation active:', active)
        setConversationActive(active)
      }
    })

    controllerRef.current = controller
    
    // Make controller globally accessible for debugging
    ;(window as any).voiceController = controller

    // Start the system
    controller.startConversation().catch((err) => {
      console.error('VoiceSystem: ❌ Failed to start conversation:', err)
      setError(`Failed to start: ${err.message}`)
    })

    return () => {
      console.log('VoiceSystem: 🧹 Cleaning up controller...')
      if (controllerRef.current) {
        controllerRef.current.emergencyReset()
        controllerRef.current = null
      }
      delete (window as any).voiceController
    }
  }, [])

  // Poll system status for real-time updates
  useEffect(() => {
    const interval = setInterval(() => {
      if (controllerRef.current) {
        const status = controllerRef.current.getSystemStatus()
        setSystemStatus(status)
        
        // Sync error state
        const currentError = controllerRef.current.getError()
        if (currentError !== error) {
          setError(currentError)
        }
      }
    }, 100)

    return () => clearInterval(interval)
  }, [error])

  // Actions
  const startConversation = useCallback(async () => {
    console.log('VoiceSystem: 🚀 Starting conversation...')
    if (controllerRef.current) {
      try {
        await controllerRef.current.startConversation()
        setError(null) // Clear any previous errors
      } catch (err: any) {
        console.error('VoiceSystem: ❌ Start conversation failed:', err)
        setError(`Start failed: ${err.message}`)
      }
    }
  }, [])

  const handleCommand = useCallback(async (command: string) => {
    console.log('VoiceSystem: 📝 Processing command:', command)
    if (controllerRef.current) {
      try {
        await controllerRef.current.handleCommand(command)
        setTranscript('') // Clear transcript after processing
      } catch (err: any) {
        console.error('VoiceSystem: ❌ Command handling failed:', err)
        setError(`Command failed: ${err.message}`)
      }
    }
  }, [])

  const handleTTSResponse = useCallback(async (text: string) => {
    console.log('VoiceSystem: 🔊 Handling TTS response:', text.substring(0, 50))
    if (controllerRef.current) {
      try {
        await controllerRef.current.handleTTSResponse(text)
      } catch (err: any) {
        console.error('VoiceSystem: ❌ TTS handling failed:', err)
        setError(`TTS failed: ${err.message}`)
      }
    }
  }, [])

  const emergencyReset = useCallback(async () => {
    console.log('VoiceSystem: 🚨 Emergency reset triggered...')
    if (controllerRef.current) {
      try {
        await controllerRef.current.emergencyReset()
        setError(null)
        setTranscript('')
        console.log('VoiceSystem: ✅ Emergency reset completed')
      } catch (err: any) {
        console.error('VoiceSystem: ❌ Emergency reset failed:', err)
        setError(`Reset failed: ${err.message}`)
      }
    }
  }, [])

  // Manual overrides for fallback UI
  const forceStartListening = useCallback(async () => {
    console.log('VoiceSystem: 🎯 Force start listening (manual override)...')
    if (controllerRef.current) {
      try {
        // Use the emergency reset followed by start conversation
        await controllerRef.current.emergencyReset()
        await new Promise(resolve => setTimeout(resolve, 500)) // Brief pause
        await controllerRef.current.startConversation()
      } catch (err: any) {
        console.error('VoiceSystem: ❌ Force start failed:', err)
        setError(`Force start failed: ${err.message}`)
      }
    }
  }, [])

  const forceStop = useCallback(() => {
    console.log('VoiceSystem: 🛑 Force stop (manual override)...')
    if (controllerRef.current) {
      // Synchronous stop - don't wait for promises
      controllerRef.current.emergencyReset().catch(console.error)
    }
  }, [])

  // Derived state
  const isListening = systemStatus.recognitionActive && 
    (systemState === SystemState.WAITING_FOR_ACTIVATION || systemState === SystemState.ACTIVE_LISTENING)
  
  const isSpeaking = systemStatus.ttsGenerating
  const isPlaying = systemStatus.ttsPlaying

  return {
    // State
    systemState,
    isListening,
    isSpeaking,
    isPlaying,
    transcript,
    conversationActive,
    error,
    systemStatus,
    
    // Actions
    startConversation,
    handleCommand,
    handleTTSResponse,
    emergencyReset,
    
    // Manual overrides
    forceStartListening,
    forceStop
  }
} 