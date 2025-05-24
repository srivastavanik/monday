import { useState, useRef, useCallback } from 'react'

interface TTSConfig {
  voiceId?: string
  apiKey?: string
  outputFormat?: string
  stability?: number
  similarityBoost?: number
  speed?: number
}

interface UseTextToSpeechReturn {
  isPlaying: boolean
  isSpeaking: boolean
  speak: (text: string) => Promise<void>
  stop: () => void
  error: string | null
  initializeAudio: () => Promise<boolean>
}

export const useTextToSpeech = (config: TTSConfig = {}): UseTextToSpeechReturn => {
  const [isPlaying, setIsPlaying] = useState(false)
  const [isSpeaking, setIsSpeaking] = useState(false)
  const [error, setError] = useState<string | null>(null)
  
  const wsRef = useRef<WebSocket | null>(null)
  const audioContextRef = useRef<AudioContext | null>(null)
  const audioQueueRef = useRef<ArrayBuffer[]>([])
  const isPlayingAudioRef = useRef(false)
  const audioInitializedRef = useRef(false)

  // Default configuration
  const {
    voiceId = (import.meta as any).env?.VITE_ELEVENLABS_VOICE_ID || 'XrExE9yKIg1WjnnlVkGX',
    apiKey = (import.meta as any).env?.VITE_ELEVENLABS_API_KEY || '',
    outputFormat = 'mp3_44100_128',
    stability = 0.5,
    similarityBoost = 0.8,
    speed = 1.0
  } = config

  // Initialize audio context with user gesture
  const initializeAudio = useCallback(async (): Promise<boolean> => {
    if (audioInitializedRef.current && audioContextRef.current?.state === 'running') {
      return true
    }

    try {
      // Create AudioContext only after user gesture
      if (!audioContextRef.current) {
        audioContextRef.current = new (window.AudioContext || (window as any).webkitAudioContext)()
      }
      
      // Resume if suspended
      if (audioContextRef.current.state === 'suspended') {
        await audioContextRef.current.resume()
      }
      
      if (audioContextRef.current.state === 'running') {
        audioInitializedRef.current = true
        console.log('AudioContext initialized successfully:', audioContextRef.current.state)
        setError(null)
        return true
      } else {
        throw new Error(`AudioContext state: ${audioContextRef.current.state}`)
      }
    } catch (err) {
      console.error('Failed to initialize audio context:', err)
      setError('Audio initialization failed - please interact with the page first')
      audioInitializedRef.current = false
      return false
    }
  }, [])

  // Play audio buffer
  const playAudioBuffer = useCallback(async (audioBuffer: ArrayBuffer): Promise<void> => {
    if (!audioContextRef.current || audioContextRef.current.state !== 'running') {
      throw new Error('AudioContext not ready')
    }

    try {
      const audioData = await audioContextRef.current.decodeAudioData(audioBuffer.slice(0))
      const source = audioContextRef.current.createBufferSource()
      source.buffer = audioData
      source.connect(audioContextRef.current.destination)
      
      return new Promise<void>((resolve) => {
        source.onended = () => resolve()
        source.start()
      })
    } catch (err) {
      console.error('Failed to play audio:', err)
      throw err
    }
  }, [])

  // Process audio queue
  const processAudioQueue = useCallback(async () => {
    if (isPlayingAudioRef.current || audioQueueRef.current.length === 0) return

    isPlayingAudioRef.current = true
    setIsPlaying(true)

    try {
      while (audioQueueRef.current.length > 0) {
        const audioBuffer = audioQueueRef.current.shift()!
        await playAudioBuffer(audioBuffer)
      }
    } catch (err) {
      console.error('Error processing audio queue:', err)
      setError('Audio playback failed')
    } finally {
      isPlayingAudioRef.current = false
      setIsPlaying(false)
    }
  }, [playAudioBuffer])

  // Convert base64 to array buffer
  const base64ToArrayBuffer = useCallback((base64: string): ArrayBuffer => {
    const binaryString = window.atob(base64)
    const len = binaryString.length
    const bytes = new Uint8Array(len)
    for (let i = 0; i < len; i++) {
      bytes[i] = binaryString.charCodeAt(i)
    }
    return bytes.buffer
  }, [])

  // Main speak function
  const speak = useCallback(async (text: string): Promise<void> => {
    if (!apiKey) {
      setError('ElevenLabs API key not configured')
      return
    }

    if (!text.trim()) return

    try {
      setError(null)
      setIsSpeaking(true)
      
      // Try to initialize audio - if it fails, we'll continue without audio
      const audioReady = await initializeAudio()
      if (!audioReady) {
        console.warn('Audio not ready, TTS will generate but not play')
        // Continue with TTS generation but without audio playback
      }

      // Create WebSocket connection
      const wsUrl = `wss://api.elevenlabs.io/v1/text-to-speech/${voiceId}/stream-input?output_format=${outputFormat}&auto_mode=true`
      const ws = new WebSocket(wsUrl)
      wsRef.current = ws

      return new Promise<void>((resolve, reject) => {
        let hasFinished = false

        ws.onopen = () => {
          console.log('ElevenLabs WebSocket connected')
          
          // Send initial configuration
          ws.send(JSON.stringify({
            text: ' ',
            voice_settings: {
              stability,
              similarity_boost: similarityBoost,
              speed
            },
            xi_api_key: apiKey
          }))

          // Send the actual text
          ws.send(JSON.stringify({
            text: text + ' ',
            try_trigger_generation: false
          }))

          // End the stream
          ws.send(JSON.stringify({
            text: ''
          }))
        }

        ws.onmessage = (event) => {
          try {
            const data = JSON.parse(event.data)
            
            if (data.audio && audioReady) {
              // Only process audio if AudioContext is ready
              const audioBuffer = base64ToArrayBuffer(data.audio)
              audioQueueRef.current.push(audioBuffer)
              processAudioQueue()
            }
            
            if (data.isFinal) {
              console.log('TTS generation complete')
              hasFinished = true
              setIsSpeaking(false)
              ws.close()
              resolve()
            }
          } catch (err) {
            console.error('Error processing TTS message:', err)
          }
        }

        ws.onerror = (error) => {
          console.error('ElevenLabs WebSocket error:', error)
          setError('Text-to-speech connection failed')
          setIsSpeaking(false)
          if (!hasFinished) reject(error)
        }

        ws.onclose = () => {
          console.log('ElevenLabs WebSocket closed')
          setIsSpeaking(false)
          wsRef.current = null
          if (!hasFinished) resolve() // Resolve even if no audio played
        }

        // Timeout after 30 seconds
        setTimeout(() => {
          if (!hasFinished) {
            ws.close()
            setIsSpeaking(false)
            setError('TTS timeout')
            reject(new Error('TTS timeout'))
          }
        }, 30000)
      })

    } catch (err) {
      console.error('TTS error:', err)
      setError('Text-to-speech failed')
      setIsSpeaking(false)
      throw err
    }
  }, [apiKey, voiceId, outputFormat, stability, similarityBoost, speed, initializeAudio, base64ToArrayBuffer, processAudioQueue])

  // Stop function
  const stop = useCallback(() => {
    if (wsRef.current) {
      wsRef.current.close()
      wsRef.current = null
    }
    
    // Clear audio queue but don't close AudioContext (keep it for reuse)
    audioQueueRef.current = []
    isPlayingAudioRef.current = false
    setIsPlaying(false)
    setIsSpeaking(false)
    setError(null)
  }, [])

  return {
    isPlaying,
    isSpeaking,
    speak,
    stop,
    error,
    initializeAudio
  }
} 