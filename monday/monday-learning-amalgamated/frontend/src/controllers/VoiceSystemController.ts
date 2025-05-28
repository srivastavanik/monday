// Type declarations for Web Speech API
declare global {
  interface Window {
    webkitSpeechRecognition: any
  }
}

interface SpeechRecognition extends EventTarget {
  continuous: boolean
  interimResults: boolean
  lang: string
  start(): void
  stop(): void
  abort(): void
  onstart: ((this: SpeechRecognition, ev: Event) => any) | null
  onresult: ((this: SpeechRecognition, ev: SpeechRecognitionEvent) => any) | null
  onerror: ((this: SpeechRecognition, ev: SpeechRecognitionErrorEvent) => any) | null
  onend: ((this: SpeechRecognition, ev: Event) => any) | null
}

interface SpeechRecognitionEvent extends Event {
  resultIndex: number
  results: SpeechRecognitionResultList
}

interface SpeechRecognitionErrorEvent extends Event {
  error: string
  message: string
}

interface SpeechRecognitionResultList {
  length: number
  item(index: number): SpeechRecognitionResult
  [index: number]: SpeechRecognitionResult
}

interface SpeechRecognitionResult {
  length: number
  item(index: number): SpeechRecognitionAlternative
  [index: number]: SpeechRecognitionAlternative
  isFinal: boolean
}

interface SpeechRecognitionAlternative {
  transcript: string
  confidence: number
}

enum SystemState {
  IDLE = 'IDLE',
  WAITING_FOR_ACTIVATION = 'WAITING_FOR_ACTIVATION',
  PROCESSING_COMMAND = 'PROCESSING_COMMAND',
  PLAYING_TTS = 'PLAYING_TTS',
  ACTIVE_LISTENING = 'ACTIVE_LISTENING',
  ERROR = 'ERROR',
  RESETTING = 'RESETTING'
}

interface BrowserFixes {
  recognitionRestartDelay: number
  abortMethod: (rec: SpeechRecognition) => void
  requiresUserInteraction: boolean
  forceSingleShot: boolean
}

interface SystemStatus {
  recognitionActive: boolean
  ttsPlaying: boolean
  ttsGenerating: boolean
  lastTransition: number
  errorCount: number
  lastError: string | null
}

interface TTSConfig {
  voiceId?: string
  apiKey?: string
  outputFormat?: string
  stability?: number
  similarityBoost?: number
  speed?: number
}

import { CommandProcessor } from './CommandProcessor'

class VoiceSystemController {
  private state: SystemState = SystemState.IDLE
  private recognition: SpeechRecognition | null = null
  private stateTransitionQueue: Array<() => Promise<void>> = []
  private transitionInProgress = false
  private audioContext: AudioContext | null = null
  private currentTranscript = ''
  
  // Event callbacks
  private onStateChange?: (state: SystemState) => void
  private onTranscriptChange?: (transcript: string) => void
  private onError?: (error: string) => void
  
  // TTS WebSocket
  private ttsWebSocket: WebSocket | null = null
  private allAudioBuffers: ArrayBuffer[] = []
  private ttsConfig: TTSConfig
  
  // Command processor for handling all command logic
  private commandProcessor = CommandProcessor.getInstance()
  
  // Single source of truth for all states
  private systemStatus: SystemStatus = {
    recognitionActive: false,
    ttsPlaying: false,
    ttsGenerating: false,
    lastTransition: Date.now(),
    errorCount: 0,
    lastError: null
  }

  constructor(ttsConfig: TTSConfig = {}) {
    this.ttsConfig = {
      voiceId: 'pNInz6obpgDQGcFmaJgB', // Adam voice
      apiKey: (import.meta as any).env?.VITE_ELEVENLABS_API_KEY || '',
      outputFormat: 'mp3_44100_128',
      stability: 0.6,
      similarityBoost: 0.8,
      speed: 1.0,
      ...ttsConfig
    }
    
    // Start the watchdog immediately
    this.startWatchdog()
    
    console.log('VoiceController: 🏗️ Initialized with TTS config:', {
      hasApiKey: !!this.ttsConfig.apiKey,
      voiceId: this.ttsConfig.voiceId,
      outputFormat: this.ttsConfig.outputFormat
    })
  }

  // ============ STATE MANAGEMENT ============
  
  public getState(): SystemState {
    return this.state
  }
  
  public getSystemStatus(): SystemStatus {
    return { ...this.systemStatus }
  }
  
  public getError(): string | null {
    return this.systemStatus.lastError
  }
  
  public getCurrentTranscript(): string {
    return this.currentTranscript
  }

  // ============ EVENT HANDLERS ============
  
  public setCallbacks(callbacks: {
    onStateChange?: (state: SystemState) => void
    onTranscriptChange?: (transcript: string) => void
    onError?: (error: string) => void
  }) {
    this.onStateChange = callbacks.onStateChange
    this.onTranscriptChange = callbacks.onTranscriptChange
    this.onError = callbacks.onError
  }

  // ============ BROWSER-SPECIFIC FIXES ============
  
  private getBrowserSpecificFixes(): BrowserFixes {
    const isChrome = /Chrome/.test(navigator.userAgent)
    const isEdge = /Edg/.test(navigator.userAgent)
    const isFirefox = /Firefox/.test(navigator.userAgent)
    
    return {
      recognitionRestartDelay: isChrome ? 500 : 300,
      abortMethod: isEdge ? 
        (rec: SpeechRecognition) => { rec.abort(); rec.abort(); } : 
        (rec: SpeechRecognition) => rec.abort(),
      requiresUserInteraction: !document.hasFocus(),
      forceSingleShot: isFirefox
    }
  }

  // ============ AUDIO CONTEXT MANAGEMENT ============
  
  private async initializeAudioContext(): Promise<boolean> {
    try {
      if (!this.audioContext || this.audioContext.state === 'closed') {
        this.audioContext = new (window.AudioContext || (window as any).webkitAudioContext)()
      }
      
      if (this.audioContext.state === 'suspended') {
        await this.audioContext.resume()
      }
      
      return this.audioContext.state === 'running'
    } catch (error) {
      console.error('VoiceController: Audio context initialization failed:', error)
      return false
    }
  }

  // ============ SPEECH RECOGNITION MANAGEMENT ============
  
  private async ensureRecognitionActive(): Promise<boolean> {
    console.log('VoiceController: 🔄 Ensuring recognition active...')
    
    // 1. ALWAYS create fresh recognition instance
    if (this.recognition) {
      const fixes = this.getBrowserSpecificFixes()
      try {
        fixes.abortMethod(this.recognition)
      } catch (e) {
        console.log('VoiceController: Ignoring cleanup error:', e)
      }
      this.recognition = null
    }
    
    // 2. Wait for browser to release resources
    const fixes = this.getBrowserSpecificFixes()
    await new Promise(resolve => setTimeout(resolve, fixes.recognitionRestartDelay))
    
    // 3. Create new instance with fresh config
    try {
      this.recognition = new (window as any).webkitSpeechRecognition()
      if (!this.recognition) {
        throw new Error('Failed to create recognition instance')
      }
      
      this.recognition.continuous = !fixes.forceSingleShot
      this.recognition.interimResults = true
      this.recognition.lang = 'en-US'
      
      // 4. Implement bulletproof event handlers
      let startSuccess = false
      
      this.recognition.onstart = () => {
        console.log('VoiceController: ✅ Recognition started successfully')
        startSuccess = true
        this.systemStatus.recognitionActive = true
        this.updateLastTransition()
      }
      
      this.recognition.onresult = (event: SpeechRecognitionEvent) => {
        let finalTranscript = ''
        for (let i = event.resultIndex; i < event.results.length; i++) {
          if (event.results[i].isFinal) {
            finalTranscript += event.results[i][0].transcript
          }
        }
        
        if (finalTranscript) {
          console.log('VoiceController: 🎤 Final transcript:', finalTranscript)
          this.currentTranscript = finalTranscript
          this.onTranscriptChange?.(finalTranscript)
          
          // Send directly to command processor - no filtering here
          this.commandProcessor.queueCommand(finalTranscript, Date.now())
        }
      }
      
      this.recognition.onerror = (event: SpeechRecognitionErrorEvent) => {
        console.error('VoiceController: ❌ Recognition error:', event.error)
        this.systemStatus.recognitionActive = false
        this.systemStatus.errorCount++
        this.systemStatus.lastError = `Recognition error: ${event.error}`
        this.onError?.(this.systemStatus.lastError)
      }
      
      this.recognition.onend = () => {
        console.log('VoiceController: 🛑 Recognition ended')
        this.systemStatus.recognitionActive = false
        this.updateLastTransition()
      }
      
      // 5. Try to start with timeout
      return new Promise<boolean>((resolve) => {
        const timeout = setTimeout(() => {
          if (!startSuccess) {
            console.warn('VoiceController: ⚠️ Recognition start timeout')
            resolve(false)
          }
        }, 3000) // Increased timeout
        
        try {
          if (this.recognition) {
            this.recognition.start()
          }
          
          // 6. Verify it actually started
          setTimeout(() => {
            clearTimeout(timeout)
            resolve(startSuccess)
          }, 1000) // Increased verification delay
        } catch (e) {
          console.error('VoiceController: ❌ Recognition start failed:', e)
          clearTimeout(timeout)
          resolve(false)
        }
      })
    } catch (error) {
      console.error('VoiceController: ❌ Recognition creation failed:', error)
      return false
    }
  }

  // ============ TTS MANAGEMENT ============
  
  private async playTTSWithGuaranteedCompletion(text: string): Promise<void> {
    console.log('VoiceController: 🔊 Starting TTS with guaranteed completion:', text.substring(0, 50))
    
    if (!this.ttsConfig.apiKey) {
      throw new Error('ElevenLabs API key not configured')
    }
    
    this.systemStatus.ttsGenerating = true
    this.allAudioBuffers = []
    
    try {
      // Initialize audio context
      const audioReady = await this.initializeAudioContext()
      if (!audioReady) {
        throw new Error('Audio context not ready')
      }
      
      // Generate TTS via WebSocket
      await this.generateTTSAudio(text)
      
      if (this.allAudioBuffers.length === 0) {
        console.warn('VoiceController: ⚠️ No audio buffers generated')
        return
      }
      
      // Calculate exact duration and play
      const totalDuration = await this.playAudioBuffersWithDuration(this.allAudioBuffers)
      console.log(`VoiceController: 🎵 TTS playback completed. Duration: ${totalDuration}ms`)
      
    } finally {
      this.systemStatus.ttsGenerating = false
      this.systemStatus.ttsPlaying = false
      this.allAudioBuffers = []
      if (this.ttsWebSocket) {
        this.ttsWebSocket.close()
        this.ttsWebSocket = null
      }
    }
  }
  
  private async generateTTSAudio(text: string): Promise<void> {
    return new Promise((resolve, reject) => {
      const wsUrl = `wss://api.elevenlabs.io/v1/text-to-speech/${this.ttsConfig.voiceId}/stream-input?output_format=${this.ttsConfig.outputFormat}&auto_mode=true`
      
      console.log('VoiceController: 🔗 Connecting to ElevenLabs WebSocket...')
      this.ttsWebSocket = new WebSocket(wsUrl)
      let hasFinished = false
      let connectionTimeout: number | null = null
      
      // Set connection timeout
      connectionTimeout = setTimeout(() => {
        if (!hasFinished && this.ttsWebSocket) {
          console.error('VoiceController: ❌ TTS WebSocket connection timeout')
          this.ttsWebSocket.close()
          reject(new Error('TTS connection timeout'))
        }
      }, 10000) // 10 second timeout
      
      this.ttsWebSocket.onopen = () => {
        console.log('VoiceController: 🔗 TTS WebSocket connected successfully')
        if (connectionTimeout) {
          clearTimeout(connectionTimeout)
          connectionTimeout = null
        }
        
        try {
          this.ttsWebSocket!.send(JSON.stringify({
            text: ' ',
            voice_settings: { 
              stability: this.ttsConfig.stability, 
              similarity_boost: this.ttsConfig.similarityBoost, 
              speed: this.ttsConfig.speed 
            },
            xi_api_key: this.ttsConfig.apiKey
          }))
          this.ttsWebSocket!.send(JSON.stringify({ text: text + ' ', try_trigger_generation: true }))
          this.ttsWebSocket!.send(JSON.stringify({ text: '' }))
          console.log('VoiceController: 📤 TTS generation request sent')
        } catch (sendError) {
          console.error('VoiceController: ❌ Failed to send TTS request:', sendError)
          reject(new Error('Failed to send TTS request'))
        }
      }
      
      this.ttsWebSocket.onmessage = (event) => {
        try {
          const data = JSON.parse(event.data)
          
          if (data.audio) {
            const audioBuffer = this.base64ToArrayBuffer(data.audio)
            this.allAudioBuffers.push(audioBuffer)
            console.log(`VoiceController: 🎵 Audio chunk ${this.allAudioBuffers.length} received`)
          }
          
          if (data.isFinal === true) {
            console.log(`VoiceController: ✅ TTS generation complete. ${this.allAudioBuffers.length} chunks received`)
            hasFinished = true
            if (connectionTimeout) {
              clearTimeout(connectionTimeout)
              connectionTimeout = null
            }
            resolve()
          }
          
          if (data.error) {
            console.error('VoiceController: ❌ ElevenLabs API error:', data.error)
            hasFinished = true
            if (connectionTimeout) {
              clearTimeout(connectionTimeout)
              connectionTimeout = null
            }
            reject(new Error(`ElevenLabs API Error: ${data.error}`))
          }
        } catch (err) {
          console.error('VoiceController: ❌ TTS message processing error:', err)
        }
      }
      
      this.ttsWebSocket.onerror = (event) => {
        console.error('VoiceController: ❌ TTS WebSocket error:', event)
        if (connectionTimeout) {
          clearTimeout(connectionTimeout)
          connectionTimeout = null
        }
        reject(new Error('TTS WebSocket connection error'))
      }
      
      this.ttsWebSocket.onclose = (event) => {
        console.log(`VoiceController: 🔌 TTS WebSocket closed. Code: ${event.code}, Reason: "${event.reason}"`)
        if (connectionTimeout) {
          clearTimeout(connectionTimeout)
          connectionTimeout = null
        }
        
        if (!hasFinished) {
          // Provide specific error messages based on close codes
          if (event.code === 1008) {
            reject(new Error(`ElevenLabs API Error: ${event.reason || 'Invalid request'}`))
          } else if (event.code === 1005) {
            reject(new Error('ElevenLabs connection failed - please check API key and internet connection'))
          } else if (event.code === 1006) {
            reject(new Error('ElevenLabs connection lost unexpectedly'))
          } else if (event.reason && event.reason.includes('quota')) {
            reject(new Error(`ElevenLabs API quota exceeded: ${event.reason}`))
          } else if (event.reason && event.reason.includes('auth')) {
            reject(new Error(`ElevenLabs authentication failed: ${event.reason}`))
          } else {
            reject(new Error(`TTS WebSocket closed prematurely (Code: ${event.code})`))
          }
        }
      }
    })
  }
  
  private async playAudioBuffersWithDuration(audioBuffers: ArrayBuffer[]): Promise<number> {
    if (!this.audioContext) {
      throw new Error('Audio context not available')
    }
    
    this.systemStatus.ttsPlaying = true
    const startTime = Date.now()
    let totalDuration = 0
    
    // Calculate total duration first
    for (const buffer of audioBuffers) {
      try {
        const audioBuffer = await this.audioContext.decodeAudioData(buffer.slice(0))
        totalDuration += audioBuffer.duration * 1000 // Convert to milliseconds
      } catch (error) {
        console.error('VoiceController: ❌ Audio decode error:', error)
      }
    }
    
    // Play all buffers sequentially
    for (const buffer of audioBuffers) {
      try {
        const audioBuffer = await this.audioContext.decodeAudioData(buffer.slice(0))
        const source = this.audioContext.createBufferSource()
        source.buffer = audioBuffer
        source.connect(this.audioContext.destination)
        
        await new Promise<void>((resolve) => {
          source.onended = () => resolve()
          source.start()
        })
      } catch (error) {
        console.error('VoiceController: ❌ Audio playback error:', error)
      }
    }
    
    // Ensure minimum duration has passed
    const elapsedTime = Date.now() - startTime
    const remainingTime = Math.max(0, totalDuration - elapsedTime + 200) // 200ms buffer
    
    if (remainingTime > 0) {
      await new Promise(resolve => setTimeout(resolve, remainingTime))
    }
    
    return totalDuration
  }
  
  private base64ToArrayBuffer(base64: string): ArrayBuffer {
    const binaryString = window.atob(base64)
    const len = binaryString.length
    const bytes = new Uint8Array(len)
    for (let i = 0; i < len; i++) {
      bytes[i] = binaryString.charCodeAt(i)
    }
    return bytes.buffer
  }

  // ============ STATE TRANSITION QUEUE ============
  
  private async transitionTo(newState: SystemState, action: () => Promise<void>): Promise<void> {
    return new Promise((resolve, reject) => {
      this.stateTransitionQueue.push(async () => {
        try {
          console.log(`VoiceController: 🔄 Transitioning ${this.state} → ${newState}`)
          const oldState = this.state
          this.state = SystemState.RESETTING
          this.onStateChange?.(this.state)
          
          await action()
          
          this.state = newState
          this.updateLastTransition()
          this.onStateChange?.(this.state)
          
          console.log(`VoiceController: ✅ Transition complete: ${oldState} → ${newState}`)
          resolve()
        } catch (e) {
          console.error(`VoiceController: ❌ Transition failed: ${this.state} → ${newState}:`, e)
          this.state = SystemState.ERROR
          this.systemStatus.lastError = `Transition failed: ${e}`
          this.systemStatus.errorCount++
          this.onStateChange?.(this.state)
          this.onError?.(this.systemStatus.lastError)
          reject(e)
        }
      })
      
      this.processTransitionQueue()
    })
  }
  
  private async processTransitionQueue(): Promise<void> {
    if (this.transitionInProgress || this.stateTransitionQueue.length === 0) {
      return
    }
    
    this.transitionInProgress = true
    const transition = this.stateTransitionQueue.shift()!
    
    try {
      await transition()
    } finally {
      this.transitionInProgress = false
      // Process next in queue
      if (this.stateTransitionQueue.length > 0) {
        setTimeout(() => this.processTransitionQueue(), 10)
      }
    }
  }

  // ============ PUBLIC API ============
  
  public async startConversation(): Promise<void> {
    console.log('VoiceController: 🚀 Starting conversation...')
    
    await this.transitionTo(SystemState.WAITING_FOR_ACTIVATION, async () => {
      const success = await this.ensureRecognitionActive()
      if (!success) {
        throw new Error('Failed to activate recognition')
      }
    })
  }
  
  public async handleCommand(command: string): Promise<void> {
    console.log('VoiceController: 📝 Handling command:', command)
    this.currentTranscript = ''
    
    await this.transitionTo(SystemState.PROCESSING_COMMAND, async () => {
      // Stop recognition during command processing
      if (this.recognition) {
        const fixes = this.getBrowserSpecificFixes()
        fixes.abortMethod(this.recognition)
        this.recognition = null
        this.systemStatus.recognitionActive = false
      }
    })
  }
  
  public async handleTTSResponse(text: string): Promise<void> {
    console.log('VoiceController: 🔊 Handling TTS response...')
    
    try {
      await this.transitionTo(SystemState.PLAYING_TTS, async () => {
        await this.playTTSWithGuaranteedCompletion(text)
      })
      
      // TTS succeeded, transition to active listening
      await this.transitionTo(SystemState.ACTIVE_LISTENING, async () => {
        // Extra delay to ensure audio device is released
        await new Promise(resolve => setTimeout(resolve, 500))
        
        // Try up to 3 times to restart recognition
        let attempts = 0
        let success = false
        
        while (attempts < 3 && !success) {
          success = await this.ensureRecognitionActive()
          if (!success) {
            attempts++
            console.warn(`VoiceController: ⚠️ Recognition restart attempt ${attempts} failed`)
            await new Promise(resolve => setTimeout(resolve, 500 * attempts))
          }
        }
        
        if (!success) {
          console.warn('VoiceController: ⚠️ Failed to restart recognition, but continuing conversation')
          // Don't throw error - just log warning and continue
        }
      })
      
    } catch (ttsError) {
      console.warn('VoiceController: ⚠️ TTS failed, continuing conversation without audio:', ttsError)
      
      // TTS failed, but maintain conversation state and continue listening
      // This is critical - don't reset conversation on TTS failure
      try {
        await this.transitionTo(SystemState.ACTIVE_LISTENING, async () => {
          // Extra delay to ensure audio device is released
          await new Promise(resolve => setTimeout(resolve, 800))
          
          // Try to restart recognition despite TTS failure
          let attempts = 0
          let success = false
          
          while (attempts < 3 && !success) {
            success = await this.ensureRecognitionActive()
            if (!success) {
              attempts++
              console.warn(`VoiceController: ⚠️ Recognition restart after TTS failure, attempt ${attempts}`)
              await new Promise(resolve => setTimeout(resolve, 700 * attempts))
            }
          }
          
          if (!success) {
            console.warn('VoiceController: ⚠️ Recognition restart failed, but keeping conversation active')
            // Even if recognition fails, keep conversation active so user can try manual restart
          }
        })
      } catch (recoveryError) {
        console.error('VoiceController: ❌ Recovery after TTS failure also failed:', recoveryError)
        // Force transition to WAITING_FOR_ACTIVATION instead of ERROR to keep system usable
        this.state = SystemState.WAITING_FOR_ACTIVATION
        this.updateLastTransition()
        this.onStateChange?.(this.state)
      }
      
      // Show user-friendly error but keep conversation going
      this.systemStatus.lastError = 'TTS temporarily unavailable, but conversation continues'
      this.onError?.('TTS temporarily unavailable, but conversation continues')
    }
  }

  // ============ EMERGENCY RECOVERY ============
  
  public async emergencyReset(): Promise<void> {
    console.log('VoiceController: 🚨 EMERGENCY RECOVERY INITIATED')
    
    // Preserve conversation state if it was active
    const wasInConversation = this.commandProcessor.isConversationActive()
    
    // Clear transition queue and reset flags
    this.stateTransitionQueue = []
    this.transitionInProgress = false
    
    // Kill everything gracefully
    if (this.recognition) {
      try { 
        const fixes = this.getBrowserSpecificFixes()
        fixes.abortMethod(this.recognition)
      } catch (e) {}
      this.recognition = null
    }
    
    if (this.ttsWebSocket) {
      try { this.ttsWebSocket.close() } catch (e) {}
      this.ttsWebSocket = null
    }
    
    if (this.audioContext) {
      try { await this.audioContext.close() } catch (e) {}
      this.audioContext = null
    }
    
    // Reset speech synthesis
    // Removed: window.speechSynthesis.cancel() - only using ElevenLabs
    
    // Clear system status but preserve conversation state
    this.systemStatus = {
      recognitionActive: false,
      ttsPlaying: false,
      ttsGenerating: false,
      lastTransition: Date.now(),
      errorCount: 0,
      lastError: null
    }
    
    this.currentTranscript = ''
    this.allAudioBuffers = []
    
    // Wait for browser to recover
    await new Promise(resolve => setTimeout(resolve, 1000))
    
    // Request fresh microphone permission
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ audio: true })
      stream.getTracks().forEach(track => track.stop())
      console.log('VoiceController: ✅ Microphone permission confirmed')
    } catch (e) {
      console.error('VoiceController: ❌ Microphone permission lost:', e)
      this.systemStatus.lastError = 'Microphone permission lost'
      // If microphone permission is lost, we have to reset conversation
      this.commandProcessor.setConversationActive(false)
    }
    
    // Restart from clean state
    this.state = SystemState.IDLE
    this.onStateChange?.(this.state)
    
    // Restore conversation state if it was active and we have microphone access
    if (wasInConversation && !this.systemStatus.lastError) {
      this.commandProcessor.setConversationActive(true)
      console.log('VoiceController: 💬 Conversation state preserved after emergency recovery')
      
      // Automatically restart conversation
      try {
        await this.startConversation()
        console.log('VoiceController: 🔄 Conversation automatically restarted after recovery')
      } catch (restartError) {
        console.warn('VoiceController: ⚠️ Failed to auto-restart conversation after recovery:', restartError)
        // Set to waiting state so user can manually restart
        this.state = SystemState.WAITING_FOR_ACTIVATION
        this.onStateChange?.(this.state)
      }
    } else {
      this.commandProcessor.setConversationActive(false)
      // Start in waiting state for new conversation
      try {
        await this.startConversation()
      } catch (startError) {
        console.warn('VoiceController: ⚠️ Failed to start conversation after recovery:', startError)
      }
    }
    
    console.log('VoiceController: ✅ Emergency recovery completed')
  }

  // ============ DIAGNOSTICS ============
  
  public detectStuckState(): string | null {
    const timeSinceLastTransition = Date.now() - this.systemStatus.lastTransition
    
    // Detect various stuck states
    if (this.state === SystemState.RESETTING && timeSinceLastTransition > 10000) {
      return 'Stuck in RESETTING state for >10s'
    }
    
    if (this.state === SystemState.PLAYING_TTS && timeSinceLastTransition > 30000) {
      return 'Stuck in PLAYING_TTS state for >30s'
    }
    
    if (this.state === SystemState.PROCESSING_COMMAND && timeSinceLastTransition > 15000) {
      return 'Stuck in PROCESSING_COMMAND state for >15s'
    }
    
    if (this.systemStatus.errorCount > 5) {
      return `Too many errors: ${this.systemStatus.errorCount}`
    }
    
    if (this.stateTransitionQueue.length > 10) {
      return `Transition queue backed up: ${this.stateTransitionQueue.length} items`
    }
    
    return null
  }
  
  private updateLastTransition(): void {
    this.systemStatus.lastTransition = Date.now()
  }
  
  private startWatchdog(): void {
    setInterval(() => {
      const stuck = this.detectStuckState()
      if (stuck) {
        console.error('VoiceController: 🚨 STUCK STATE DETECTED:', stuck)
        // Auto-recovery for certain conditions
        if (stuck.includes('Stuck in') || stuck.includes('Too many errors')) {
          console.log('VoiceController: 🔧 Attempting auto-recovery...')
          this.emergencyReset().catch(console.error)
        }
      }
    }, 5000)
  }
}

export { VoiceSystemController, SystemState }
export type { SystemStatus, TTSConfig } 