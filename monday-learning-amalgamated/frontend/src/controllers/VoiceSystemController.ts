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

class VoiceSystemController {
  private state: SystemState = SystemState.IDLE
  private recognition: SpeechRecognition | null = null
  private stateTransitionQueue: Array<() => Promise<void>> = []
  private transitionInProgress = false
  private audioContext: AudioContext | null = null
  private currentTranscript = ''
  private conversationActive = false
  
  // Event callbacks
  private onStateChange?: (state: SystemState) => void
  private onTranscriptChange?: (transcript: string) => void
  private onError?: (error: string) => void
  private onConversationChange?: (active: boolean) => void
  
  // TTS WebSocket
  private ttsWebSocket: WebSocket | null = null
  private allAudioBuffers: ArrayBuffer[] = []
  private ttsConfig: TTSConfig
  
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
      voiceId: (import.meta as any).env?.VITE_ELEVENLABS_VOICE_ID || 'XrExE9yKIg1WjnnlVkGX',
      apiKey: (import.meta as any).env?.VITE_ELEVENLABS_API_KEY || '',
      outputFormat: 'mp3_44100_128',
      stability: 0.6,
      similarityBoost: 0.8,
      speed: 1.1,
      ...ttsConfig
    }
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
  
  public isConversationActive(): boolean {
    return this.conversationActive
  }

  // ============ EVENT HANDLERS ============
  
  public setCallbacks(callbacks: {
    onStateChange?: (state: SystemState) => void
    onTranscriptChange?: (transcript: string) => void
    onError?: (error: string) => void
    onConversationChange?: (active: boolean) => void
  }) {
    this.onStateChange = callbacks.onStateChange
    this.onTranscriptChange = callbacks.onTranscriptChange
    this.onError = callbacks.onError
    this.onConversationChange = callbacks.onConversationChange
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
  
  /**
   * Stops any currently playing TTS audio and resets relevant state.
   */
  private async stopTTS(): Promise<void> {
    console.log('VoiceController: 🛑 Stopping TTS and cleaning up...')
    
    // Stop recognition immediately to prevent feedback
    if (this.recognition) {
      try {
        const fixes = this.getBrowserSpecificFixes()
        fixes.abortMethod(this.recognition)
        this.systemStatus.recognitionActive = false
        console.log('VoiceController: 🎤 Recognition stopped for TTS cleanup')
      } catch (e) {
        console.log('VoiceController: Ignoring recognition cleanup error:', e)
      }
    }
    
    // Stop and close the audio context if it exists
    if (this.audioContext) {
      try {
        await this.audioContext.close();
      } catch (e) {
        // Ignore errors
      }
      this.audioContext = null;
    }
    // Close any open TTS WebSocket
    if (this.ttsWebSocket) {
      try {
        this.ttsWebSocket.close();
      } catch (e) {}
      this.ttsWebSocket = null;
    }
    // Reset playback flags and buffers
    this.systemStatus.ttsPlaying = false;
    this.systemStatus.ttsGenerating = false;
    this.allAudioBuffers = [];
  }

  private async playTTSWithGuaranteedCompletion(text: string): Promise<void> {
    // Stop any currently playing TTS before starting new one
    await this.stopTTS();
    console.log('VoiceController: 🔊 Starting TTS with guaranteed completion:', text.substring(0, 50))
    
    if (!this.ttsConfig.apiKey) {
      console.warn('VoiceController: ⚠️ ElevenLabs API key not configured, skipping TTS')
      return
    }
    
    // Clear any existing audio buffers
    this.allAudioBuffers = []
    this.systemStatus.ttsGenerating = true
    this.systemStatus.ttsPlaying = false
    
    try {
      // Force create a new audio context to avoid caching issues
      if (this.audioContext) {
        try {
          await this.audioContext.close()
        } catch (e) {
          // Ignore errors
        }
        this.audioContext = null
      }
      
      // Initialize fresh audio context
      const audioReady = await this.initializeAudioContext()
      if (!audioReady) {
        throw new Error('Audio context not ready')
      }
      
      // Generate TTS via WebSocket with retry logic
      await this.generateTTSAudioWithRetry(text)
      
      if (this.allAudioBuffers.length === 0) {
        console.warn('VoiceController: ⚠️ No audio buffers generated, skipping playback')
        return
      }
      
      // Calculate exact duration and play
      const totalDuration = await this.playAudioBuffersWithDuration(this.allAudioBuffers)
      console.log(`VoiceController: 🎵 TTS playback completed. Duration: ${totalDuration}ms`)
      
      // Add extra delay to ensure audio is fully finished
      await new Promise(resolve => setTimeout(resolve, 1000))
      
    } catch (error) {
      console.error('VoiceController: ❌ TTS playback failed:', error)
      // Don't throw the error, just log it and continue
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
  
  private async generateTTSAudioWithRetry(text: string, maxRetries: number = 2): Promise<void> {
    for (let attempt = 1; attempt <= maxRetries; attempt++) {
      try {
        await this.generateTTSAudio(text)
        return // Success, exit retry loop
      } catch (error) {
        console.error(`VoiceController: ❌ TTS attempt ${attempt} failed:`, error)
        
        if (attempt === maxRetries) {
          console.warn('VoiceController: ⚠️ All TTS attempts failed, skipping audio playback')
          throw error
        }
        
        // Wait before retry
        await new Promise(resolve => setTimeout(resolve, 1000 * attempt))
      }
    }
  }
  
  private async generateTTSAudio(text: string): Promise<void> {
    return new Promise((resolve, reject) => {
      const wsUrl = `wss://api.elevenlabs.io/v1/text-to-speech/${this.ttsConfig.voiceId}/stream-input?output_format=${this.ttsConfig.outputFormat}&auto_mode=true`
      
      this.ttsWebSocket = new WebSocket(wsUrl)
      let hasFinished = false
      let connectionTimeout: number
      
      // Set connection timeout
      connectionTimeout = setTimeout(() => {
        if (!hasFinished) {
          console.warn('VoiceController: ⚠️ TTS WebSocket connection timeout')
          if (this.ttsWebSocket) {
            this.ttsWebSocket.close()
          }
          reject(new Error('TTS WebSocket connection timeout'))
        }
      }, 10000) // 10 second timeout
      
      this.ttsWebSocket.onopen = () => {
        console.log('VoiceController: 🔗 TTS WebSocket connected')
        clearTimeout(connectionTimeout)
        
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
        } catch (sendError) {
          console.error('VoiceController: ❌ Failed to send TTS data:', sendError)
          reject(sendError)
        }
      }
      
      this.ttsWebSocket.onmessage = (event) => {
        try {
          const data = JSON.parse(event.data)
          
          if (data.audio) {
            const audioBuffer = this.base64ToArrayBuffer(data.audio)
            this.allAudioBuffers.push(audioBuffer)
          }
          
          if (data.isFinal === true) {
            console.log(`VoiceController: ✅ TTS generation complete. ${this.allAudioBuffers.length} chunks received`)
            hasFinished = true
            resolve()
          }
        } catch (err) {
          console.error('VoiceController: ❌ TTS message processing error:', err)
        }
      }
      
      this.ttsWebSocket.onerror = (event) => {
        console.error('VoiceController: ❌ TTS WebSocket error:', event)
        clearTimeout(connectionTimeout)
        if (!hasFinished) {
          reject(new Error('TTS WebSocket error'))
        }
      }
      
      this.ttsWebSocket.onclose = (event) => {
        console.log(`VoiceController: 🔌 TTS WebSocket closed. Code: ${event.code}, Reason: "${event.reason}"`)
        clearTimeout(connectionTimeout)
        
        if (!hasFinished) {
          if (event.code === 1008 || event.reason.includes('quota')) {
            reject(new Error(`ElevenLabs API Error: ${event.reason}`))
          } else if (event.code === 1005 || event.code === 1006) {
            // Connection closed abnormally - likely rate limit or network issue
            reject(new Error('TTS connection closed unexpectedly'))
          } else {
            reject(new Error('TTS WebSocket closed prematurely'))
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
    try {
      console.log(`VoiceController: 🔄 Transitioning ${this.state} → ${newState}`)
      const oldState = this.state
      this.state = newState
      this.onStateChange?.(this.state)
      
      await action()
      
      console.log(`VoiceController: ✅ Transition complete: ${oldState} → ${newState}`)
    } catch (e) {
      console.error(`VoiceController: ❌ Transition failed: ${this.state} → ${newState}:`, e)
      this.state = SystemState.ERROR
      this.systemStatus.lastError = `Transition failed: ${e}`
      this.systemStatus.errorCount++
      this.onStateChange?.(this.state)
      this.onError?.(this.systemStatus.lastError)
      throw e
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
    
    // Set conversation active after any interaction
    if (!this.conversationActive) {
      this.conversationActive = true
      this.onConversationChange?.(true)
    }
    
    try {
      // Play TTS and wait for it to complete
      await this.transitionTo(SystemState.PLAYING_TTS, async () => {
        await this.playTTSWithGuaranteedCompletion(text)
      })
      
      // Add a longer delay after TTS completes to prevent feedback
      console.log('VoiceController: ⏳ Waiting for audio to fully complete before restarting microphone...')
      await new Promise(resolve => setTimeout(resolve, 2000)) // Increased to 2 seconds
      
      // Transition to active listening
      await this.transitionTo(SystemState.ACTIVE_LISTENING, async () => {
        // Stop any existing recognition
        if (this.recognition) {
          const fixes = this.getBrowserSpecificFixes()
          fixes.abortMethod(this.recognition)
          this.recognition = null
          this.systemStatus.recognitionActive = false
        }
        
        // Wait a bit more before starting new recognition
        await new Promise(resolve => setTimeout(resolve, 500))
        
        // Create new recognition instance
        const SpeechRecognition = window.webkitSpeechRecognition || (window as any).SpeechRecognition
        if (!SpeechRecognition) {
          throw new Error('Speech recognition not supported')
        }
        
        this.recognition = new SpeechRecognition()
        if (!this.recognition) {
          throw new Error('Failed to create speech recognition instance')
        }
        
        this.recognition.continuous = true
        this.recognition.interimResults = true
        this.recognition.lang = 'en-US'
        
        // Set up event handlers
        this.recognition.onstart = () => {
          console.log('VoiceController: ✅ Recognition restarted after TTS')
          this.systemStatus.recognitionActive = true
        }
        
        this.recognition.onresult = (event: SpeechRecognitionEvent) => {
          const result = event.results[event.results.length - 1]
          if (result.isFinal) {
            this.currentTranscript = result[0].transcript
            console.log('VoiceController: 🎤 Final transcript:', this.currentTranscript)
            this.onTranscriptChange?.(this.currentTranscript)
          }
        }
        
        this.recognition.onerror = (event: SpeechRecognitionErrorEvent) => {
          console.error('VoiceController: ❌ Recognition error:', event.error)
          this.onError?.(`Recognition error: ${event.error}`)
        }
        
        this.recognition.onend = () => {
          console.log('VoiceController: 🛑 Recognition ended')
          this.systemStatus.recognitionActive = false
          
          // Restart recognition if we're still supposed to be listening
          if (this.state === SystemState.ACTIVE_LISTENING && this.recognition) {
            setTimeout(() => {
              if (this.recognition && this.state === SystemState.ACTIVE_LISTENING) {
                try {
                  this.recognition.start()
                } catch (e) {
                  console.warn('VoiceController: ⚠️ Failed to restart recognition:', e)
                }
              }
            }, 500) // Increased restart delay
          }
        }
        
        // Start recognition
        try {
          this.recognition.start()
        } catch (error) {
          console.error('VoiceController: ❌ Failed to start recognition:', error)
          throw new Error('Failed to start recognition')
        }
      })
    } catch (error) {
      console.error('VoiceController: ❌ Error in handleTTSResponse:', error)
      // Don't throw the error, just log it and continue
    }

    // Ensure conversation stays active after TTS
    if (!this.conversationActive) {
      this.conversationActive = true
      this.onConversationChange?.(true)
    }
  }

  // ============ EMERGENCY RECOVERY ============
  
  public async emergencyReset(): Promise<void> {
    console.log('VoiceController: 🚨 EMERGENCY RECOVERY INITIATED')
    
    // Kill everything
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
    if (window.speechSynthesis) {
      window.speechSynthesis.cancel()
    }
    
    // Clear all state
    this.systemStatus = {
      recognitionActive: false,
      ttsPlaying: false,
      ttsGenerating: false,
      lastTransition: Date.now(),
      errorCount: 0,
      lastError: null
    }
    
    this.currentTranscript = ''
    this.conversationActive = false
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
    }
    
    // Restart from clean state
    this.state = SystemState.IDLE
    this.onStateChange?.(this.state)
    this.onConversationChange?.(false)
    
    console.log('VoiceController: ✅ Emergency recovery completed')

    // Automatically transition back to WAITING_FOR_ACTIVATION
    try {
      await this.transitionTo(SystemState.WAITING_FOR_ACTIVATION, async () => {
        const success = await this.ensureRecognitionActive()
        if (!success) {
          throw new Error('Failed to activate recognition')
        }
      })
    } catch (e) {
      console.error('VoiceController: ❌ Failed to restart recognition after recovery:', e)
    }
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