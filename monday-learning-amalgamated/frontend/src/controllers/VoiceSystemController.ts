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
  
  // CRITICAL: Hardware-level microphone control
  private microphoneStream: MediaStream | null = null
  private isMicrophoneLocked = false
  private ttsLockTimeout: number | null = null
  private globalTTSLock = false
  private ttsLockStartTime = 0
  
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

  // ============ HARDWARE MICROPHONE CONTROL ============
  
  private async lockMicrophone(): Promise<void> {
    console.log('VoiceController: 🔒 NUCLEAR MICROPHONE SHUTDOWN')
    this.isMicrophoneLocked = true
    this.globalTTSLock = true
    this.ttsLockStartTime = Date.now()
    
    // Set a MASSIVE timeout to force unlock if something goes wrong
    if (this.ttsLockTimeout) {
      clearTimeout(this.ttsLockTimeout)
    }
    this.ttsLockTimeout = setTimeout(() => {
      console.log('VoiceController: ⚠️ FORCE UNLOCKING microphone after timeout')
      this.isMicrophoneLocked = false
      this.globalTTSLock = false
      this.ttsLockTimeout = null
    }, 60000) // 60 second emergency timeout
    
    // NUCLEAR OPTION: Destroy EVERYTHING related to audio input
    
    // 1. Kill recognition with extreme prejudice
    if (this.recognition) {
      try {
        this.recognition.abort()
        this.recognition.stop()
        this.recognition.onstart = null
        this.recognition.onresult = null
        this.recognition.onerror = null
        this.recognition.onend = null
        this.recognition = null
        console.log('VoiceController: 💀 Recognition DESTROYED')
      } catch (e) {
        console.log('VoiceController: Ignoring recognition destruction error:', e)
      }
    }
    
    // 2. Kill ALL microphone streams
    if (this.microphoneStream) {
      this.microphoneStream.getTracks().forEach(track => {
        track.stop()
        track.enabled = false
        console.log('VoiceController: 💀 Microphone track DESTROYED')
      })
      this.microphoneStream = null
    }
    
    // 3. Kill ALL audio streams from getUserMedia
    try {
      const allStreams = await navigator.mediaDevices.enumerateDevices()
      console.log('VoiceController: 💀 Attempting to kill all audio devices')
    } catch (e) {
      console.log('VoiceController: Could not enumerate devices:', e)
    }
    
    // 4. Force garbage collection if available
    if ((window as any).gc) {
      (window as any).gc()
    }
    
    this.systemStatus.recognitionActive = false
    console.log('VoiceController: ✅ NUCLEAR MICROPHONE SHUTDOWN COMPLETE - NO AUDIO INPUT POSSIBLE')
  }
  
  private async unlockMicrophone(): Promise<void> {
    console.log('VoiceController: 🔓 System restart after TTS')
    
    // Clear any existing timeout
    if (this.ttsLockTimeout) {
      clearTimeout(this.ttsLockTimeout)
      this.ttsLockTimeout = null
    }
    
    // Calculate how long the lock was active
    const lockDuration = Date.now() - this.ttsLockStartTime
    console.log(`VoiceController: 🔓 Lock was active for ${lockDuration}ms`)
    
    // Ensure reasonable minimum lock time to prevent overlap
    const minimumLockTime = 5000 // Reduced to 5 seconds
    if (lockDuration < minimumLockTime) {
      const additionalWait = minimumLockTime - lockDuration
      console.log(`VoiceController: ⏳ Additional wait of ${additionalWait}ms for audio separation`)
      await new Promise(resolve => setTimeout(resolve, additionalWait))
    }
    
    // Reasonable wait for residual audio to clear
    console.log('VoiceController: ⏳ Wait for audio system clearance...')
    await new Promise(resolve => setTimeout(resolve, 2000)) // Reduced to 2 seconds
    
    // Force browser to release audio resources
    if ((window as any).gc) {
      (window as any).gc()
      console.log('VoiceController: 🗑️ Forced garbage collection')
    }
    
    // Short delay after garbage collection
    await new Promise(resolve => setTimeout(resolve, 1000)) // Reduced to 1 second
    
    this.isMicrophoneLocked = false
    this.globalTTSLock = false
    console.log('VoiceController: ✅ System ready for restart - audio separated')
  }
  
  private async requestMicrophoneAccess(): Promise<MediaStream | null> {
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ 
        audio: {
          echoCancellation: true,
          noiseSuppression: true,
          autoGainControl: true
        } 
      })
      this.microphoneStream = stream
      console.log('VoiceController: 🎤 Fresh microphone access granted')
      return stream
    } catch (error) {
      console.error('VoiceController: ❌ Failed to get microphone access:', error)
      return null
    }
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
      
      // Check if AudioContext is actually running
      if (this.audioContext.state !== 'running') {
        console.warn('VoiceController: ⚠️ AudioContext not running, likely needs user gesture')
        return false
      }
      
      return this.audioContext.state === 'running'
    } catch (error) {
      console.warn('VoiceController: ⚠️ Audio context initialization failed (likely needs user gesture):', error)
      return false
    }
  }

  // ============ SPEECH RECOGNITION MANAGEMENT ============
  
  private async ensureRecognitionActive(): Promise<boolean> {
    console.log('VoiceController: 🔄 Ensuring recognition active...')
    
    // CRITICAL SAFETY CHECK: Never start recognition during TTS or when locked
    if (this.systemStatus.ttsPlaying || this.systemStatus.ttsGenerating || this.isMicrophoneLocked || this.globalTTSLock) {
      console.log('VoiceController: 🔇 BLOCKED recognition start - TTS active, microphone locked, or global TTS lock active')
      return false
    }
    
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
    
    // 2. Ensure we have fresh microphone access
    if (!this.microphoneStream) {
      const stream = await this.requestMicrophoneAccess()
      if (!stream) {
        console.error('VoiceController: ❌ Cannot start recognition - no microphone access')
        return false
      }
    }
    
    // 3. Wait for browser to release resources
    const fixes = this.getBrowserSpecificFixes()
    await new Promise(resolve => setTimeout(resolve, fixes.recognitionRestartDelay))
    
    // 4. Create new instance with fresh config
    try {
      this.recognition = new (window as any).webkitSpeechRecognition()
      if (!this.recognition) {
        throw new Error('Failed to create recognition instance')
      }
      
      // Improved settings for better recognition
      this.recognition.continuous = true
      this.recognition.interimResults = true
      this.recognition.lang = 'en-US'
      // maxAlternatives not supported in TypeScript interface
      
      // 5. Implement bulletproof event handlers
      let startSuccess = false
      let lastResultTime = Date.now()
      
      this.recognition.onstart = () => {
        console.log('VoiceController: ✅ Recognition started successfully')
        startSuccess = true
        this.systemStatus.recognitionActive = true
        this.updateLastTransition()
      }
      
      this.recognition.onresult = (event: SpeechRecognitionEvent) => {
        // Double-check lock status before processing any results
        if (this.isMicrophoneLocked || this.systemStatus.ttsPlaying || this.systemStatus.ttsGenerating || this.globalTTSLock) {
          console.log('VoiceController: 🔇 Ignoring recognition result - microphone locked, TTS active, or global TTS lock active')
          return
        }
        
        lastResultTime = Date.now()
        let finalTranscript = ''
        let interimTranscript = ''
        
        for (let i = event.resultIndex; i < event.results.length; i++) {
          const transcript = event.results[i][0].transcript
          if (event.results[i].isFinal) {
            finalTranscript += transcript
          } else {
            interimTranscript += transcript
          }
        }
        
        if (finalTranscript) {
          console.log('VoiceController: 🎤 Final transcript:', finalTranscript)
          this.currentTranscript = finalTranscript.trim()
          this.onTranscriptChange?.(this.currentTranscript)
        }
      }
      
      this.recognition.onerror = (event: SpeechRecognitionErrorEvent) => {
        // Only log actual errors, not "no-speech" which is normal
        if (event.error !== 'no-speech') {
          console.error('VoiceController: ❌ Recognition error:', event.error)
          this.systemStatus.recognitionActive = false
          this.systemStatus.errorCount++
          this.systemStatus.lastError = `Recognition error: ${event.error}`
          this.onError?.(this.systemStatus.lastError)
        }
      }
      
      this.recognition.onend = () => {
        console.log('VoiceController: 🛑 Recognition ended')
        this.systemStatus.recognitionActive = false
        this.updateLastTransition()
        
        // CRITICAL: Do NOT restart recognition if TTS is playing, generating, or microphone is locked
        if (this.systemStatus.ttsPlaying || this.systemStatus.ttsGenerating || this.isMicrophoneLocked || this.globalTTSLock) {
          console.log('VoiceController: 🔇 NOT restarting recognition - TTS active, microphone locked, or global TTS lock active')
          return
        }
        
        // Auto-restart if we should still be listening and it's been a while since last result
        const timeSinceLastResult = Date.now() - lastResultTime
        if (this.state === SystemState.WAITING_FOR_ACTIVATION || 
            this.state === SystemState.ACTIVE_LISTENING) {
          
          // Only restart if it's been more than 2 seconds since last speech
          if (timeSinceLastResult > 2000) {
            console.log('VoiceController: 🔄 Auto-restarting recognition after silence')
            setTimeout(() => {
              // Triple-check all conditions before restarting
              if (!this.systemStatus.ttsPlaying && !this.systemStatus.ttsGenerating && 
                  !this.isMicrophoneLocked && !this.globalTTSLock &&
                  (this.state === SystemState.WAITING_FOR_ACTIVATION || 
                   this.state === SystemState.ACTIVE_LISTENING) && 
                  !this.systemStatus.recognitionActive) {
                this.ensureRecognitionActive().catch(console.error)
              } else {
                console.log('VoiceController: 🔇 Skipping auto-restart - conditions not met (including global TTS lock)')
              }
            }, 1000)
          }
        }
      }
      
      // 6. Try to start with longer timeout
      return new Promise<boolean>((resolve) => {
        const timeout = setTimeout(() => {
          if (!startSuccess) {
            console.warn('VoiceController: ⚠️ Recognition start timeout')
            resolve(false)
          }
        }, 5000) // Increased timeout to 5 seconds
        
        try {
          if (this.recognition) {
            this.recognition.start()
          }
          
          // 7. Verify it actually started
          setTimeout(() => {
            clearTimeout(timeout)
            resolve(startSuccess)
          }, 1500) // Increased verification delay
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
    
    // IMMEDIATELY stop recognition to prevent feedback
    if (this.recognition) {
      try {
        const fixes = this.getBrowserSpecificFixes()
        fixes.abortMethod(this.recognition)
        this.recognition = null
        this.systemStatus.recognitionActive = false
        console.log('VoiceController: 🎤 Recognition STOPPED for TTS cleanup')
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
    
    // CRITICAL: Ensure microphone is completely disabled during TTS
    console.log('VoiceController: 🔇 FORCING microphone OFF for TTS')
    if (this.recognition) {
      const fixes = this.getBrowserSpecificFixes()
      try {
        fixes.abortMethod(this.recognition)
        this.recognition.abort() // Double abort for safety
        this.recognition = null
      } catch (e) {
        console.log('VoiceController: Ignoring recognition abort error:', e)
      }
    }
    this.systemStatus.recognitionActive = false
    
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
        console.warn('VoiceController: ⚠️ Audio not available (likely needs user gesture), skipping TTS playback')
        // Still wait a reasonable amount of time to simulate TTS duration
        const estimatedDuration = Math.max(2000, text.length * 50) // ~50ms per character
        await new Promise(resolve => setTimeout(resolve, estimatedDuration))
        return
      }
      
      // Generate TTS via WebSocket with retry logic
      await this.generateTTSAudioWithRetry(text)
      
      if (this.allAudioBuffers.length === 0) {
        console.warn('VoiceController: ⚠️ No audio buffers generated, skipping playback')
        return
      }
      
      // Mark TTS as playing to prevent any microphone activation
      this.systemStatus.ttsPlaying = true
      console.log('VoiceController: 🔇 TTS PLAYING - Microphone MUST stay OFF')
      
      // Calculate exact duration and play
      const totalDuration = await this.playAudioBuffersWithDuration(this.allAudioBuffers)
      console.log(`VoiceController: 🎵 TTS playback completed. Duration: ${totalDuration}ms`)
      
      // Add extra delay to ensure audio is fully finished and system audio has cleared
      console.log('VoiceController: ⏳ Post-TTS delay to ensure audio cleared...')
      await new Promise(resolve => setTimeout(resolve, 3000)) // Reduced to 3 seconds
      
    } catch (error) {
      console.warn('VoiceController: ⚠️ TTS playback failed, continuing without audio:', error)
      // Don't throw the error, just log it and continue
      // Wait a reasonable amount of time to simulate TTS duration
      const estimatedDuration = Math.max(3000, text.length * 80) // Reduced to 80ms per character
      await new Promise(resolve => setTimeout(resolve, estimatedDuration))
    } finally {
      this.systemStatus.ttsGenerating = false
      this.systemStatus.ttsPlaying = false
      this.allAudioBuffers = []
      if (this.ttsWebSocket) {
        this.ttsWebSocket.close()
        this.ttsWebSocket = null
      }
      console.log('VoiceController: ✅ TTS finished with isolation - microphone can be reactivated')
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
    console.log('VoiceController: 🔊 TTS RESPONSE WITH STRONG AUDIO ISOLATION...')
    
    // Set conversation active after any interaction
    if (!this.conversationActive) {
      this.conversationActive = true
      this.onConversationChange?.(true)
    }
    
    try {
      // STEP 1: NUCLEAR MICROPHONE SHUTDOWN
      await this.lockMicrophone()
      
      // STEP 2: Short pre-TTS delay to ensure microphone shutdown
      console.log('VoiceController: ⏳ Pre-TTS delay for microphone shutdown...')
      await new Promise(resolve => setTimeout(resolve, 1000)) // Reduced to 1 second
      
      // STEP 3: Play TTS with microphone guaranteed destroyed
      await this.transitionTo(SystemState.PLAYING_TTS, async () => {
        console.log('VoiceController: 🔇 Microphone DESTROYED - playing TTS in isolation')
        
        // Play TTS with nuclear-level protection
        await this.playTTSWithGuaranteedCompletion(text)
        
        // STEP 4: Reasonable delay to ensure audio has cleared
        console.log('VoiceController: ⏳ Post-TTS delay for audio clearance...')
        await new Promise(resolve => setTimeout(resolve, 3000)) // Reduced to 3 seconds
      })
      
      // STEP 5: Complete system restart with fresh everything
      await this.transitionTo(SystemState.ACTIVE_LISTENING, async () => {
        console.log('VoiceController: 🔄 System restart after TTS isolation')
        
        // Short pre-unlock delay
        await new Promise(resolve => setTimeout(resolve, 1000)) // Reduced to 1 second
        
        // Unlock the microphone (which includes delays)
        await this.unlockMicrophone()
        
        // Reasonable additional safety delay before recognition
        console.log('VoiceController: ⏳ Post-unlock safety delay...')
        await new Promise(resolve => setTimeout(resolve, 2000)) // Reduced to 2 seconds
        
        // Start completely fresh recognition with brand new microphone access
        const success = await this.ensureRecognitionActive()
        if (!success) {
          console.warn('VoiceController: ⚠️ Failed to restart recognition after TTS')
          this.systemStatus.lastError = 'Failed to restart recognition after TTS'
          this.systemStatus.errorCount++
        } else {
          console.log('VoiceController: ✅ Microphone reactivated with system restart')
        }
      })
      
    } catch (error) {
      console.error('VoiceController: ❌ Error in TTS response:', error)
      
      // CRITICAL: Always unlock microphone even if there's an error
      try {
        console.log('VoiceController: 🔧 Emergency microphone unlock after TTS error')
        await this.unlockMicrophone()
        await new Promise(resolve => setTimeout(resolve, 2000)) // Reduced to 2 seconds
        await this.ensureRecognitionActive()
      } catch (unlockError) {
        console.error('VoiceController: ❌ Failed to unlock microphone after TTS error:', unlockError)
        this.systemStatus.lastError = 'Critical: Failed to unlock microphone after TTS'
        this.systemStatus.errorCount++
      }
    }

    // Ensure conversation stays active after TTS
    if (!this.conversationActive) {
      this.conversationActive = true
      this.onConversationChange?.(true)
    }
  }

  public async interruptTTS(): Promise<void> {
    console.log('VoiceController: ⚡ Interrupting TTS with hardware unlock')
    
    // Stop TTS immediately and unlock microphone
    await this.stopTTS()
    await this.unlockMicrophone()
    
    // Transition to active listening
    await this.transitionTo(SystemState.ACTIVE_LISTENING, async () => {
      // Wait a moment for audio to clear
      await new Promise(resolve => setTimeout(resolve, 500))
      
      // Start fresh recognition with new microphone access
      const success = await this.ensureRecognitionActive()
      if (!success) {
        console.warn('VoiceController: ⚠️ Failed to start recognition after TTS interruption')
        throw new Error('Failed to start recognition after TTS interruption')
      }
      
      console.log('VoiceController: ✅ Microphone active after TTS interruption')
    })
  }

  // ============ EMERGENCY RECOVERY ============
  
  public async emergencyReset(): Promise<void> {
    console.log('VoiceController: 🚨 EMERGENCY RECOVERY INITIATED')
    
    // Clear any TTS lock timeout
    if (this.ttsLockTimeout) {
      clearTimeout(this.ttsLockTimeout)
      this.ttsLockTimeout = null
    }
    
    // Force unlock microphone
    this.isMicrophoneLocked = false
    this.globalTTSLock = false
    
    // Kill everything
    if (this.recognition) {
      try { 
        const fixes = this.getBrowserSpecificFixes()
        fixes.abortMethod(this.recognition)
        this.recognition.abort()
        this.recognition.stop()
      } catch (e) {}
      this.recognition = null
    }
    
    // Stop microphone stream
    if (this.microphoneStream) {
      this.microphoneStream.getTracks().forEach(track => track.stop())
      this.microphoneStream = null
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
    await new Promise(resolve => setTimeout(resolve, 2000))
    
    // Request fresh microphone permission
    try {
      const stream = await this.requestMicrophoneAccess()
      if (stream) {
        console.log('VoiceController: ✅ Fresh microphone access confirmed')
      } else {
        console.error('VoiceController: ❌ Failed to get fresh microphone access')
        this.systemStatus.lastError = 'Microphone permission lost'
      }
    } catch (e) {
      console.error('VoiceController: ❌ Microphone permission error:', e)
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