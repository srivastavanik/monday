import { useState, useEffect, useCallback, useRef } from 'react'

interface VoiceRecognitionHook {
  isListening: boolean
  transcript: string
  confidence: number
  error: string | null
  isSupported: boolean
  startListening: () => void
  stopListening: () => void
  resetTranscript: () => void
}

export const useVoiceRecognition = (): VoiceRecognitionHook => {
  const [isListening, setIsListening] = useState(false)
  const [transcript, setTranscript] = useState('')
  const [confidence, setConfidence] = useState(0)
  const [error, setError] = useState<string | null>(null)
  const [isSupported, setIsSupported] = useState(false)
  
  const recognitionRef = useRef<any>(null)
  const lastTranscriptTimeRef = useRef<number>(Date.now())

  // Check for Web Speech API support
  useEffect(() => {
    const SpeechRecognition = (window as any).SpeechRecognition || (window as any).webkitSpeechRecognition
    if (SpeechRecognition) {
      setIsSupported(true)
      recognitionRef.current = new SpeechRecognition()
    } else {
      setIsSupported(false)
      setError('Speech recognition not supported in this browser')
    }
  }, [])

  const startListening = useCallback(() => {
    if (!recognitionRef.current || !isSupported || isListening) {
      return
    }

    try {
      setError(null)
      setTranscript('')
      lastTranscriptTimeRef.current = Date.now() // Reset watchdog timer
      recognitionRef.current.start()
    } catch (err) {
      setError('Failed to start voice recognition')
    }
  }, [isListening, isSupported])

  const stopListening = useCallback(() => {
    if (recognitionRef.current && isListening) {
      recognitionRef.current.stop()
    }
  }, [isListening])

  const resetTranscript = useCallback(() => {
    setTranscript('')
    setConfidence(0)
  }, [])

  // Configure speech recognition
  useEffect(() => {
    if (!recognitionRef.current) return

    const recognition = recognitionRef.current

    recognition.continuous = true
    recognition.interimResults = true
    recognition.lang = 'en-US'

    recognition.onstart = () => {
      setIsListening(true)
      setError(null)
    }

    recognition.onresult = (event: any) => {
      let finalTranscript = ''
      let interimTranscript = ''

      for (let i = event.resultIndex; i < event.results.length; i++) {
        const result = event.results[i]
        const transcriptText = result[0].transcript

        if (result.isFinal) {
          finalTranscript += transcriptText
          setConfidence(result[0].confidence)
        } else {
          interimTranscript += transcriptText
        }
      }

      const fullTranscript = finalTranscript + interimTranscript
      setTranscript(fullTranscript)
      
      // Update last transcript time for watchdog
      if (fullTranscript.trim()) {
        lastTranscriptTimeRef.current = Date.now()
      }
    }

    recognition.onerror = (event: any) => {
      console.error('Speech recognition error:', event.error)
      setError(`Speech recognition error: ${event.error}`)
      setIsListening(false)
    }

    recognition.onend = () => {
      setIsListening(false)
    }
  }, [])

  // Watchdog timer to restart speech recognition if stuck
  useEffect(() => {
    if (!isListening || !recognitionRef.current) return

    const watchdogTimer = setInterval(() => {
      const timeSinceLastTranscript = Date.now() - lastTranscriptTimeRef.current
      
      // If listening for more than 30 seconds without any transcript, restart
      if (timeSinceLastTranscript > 30000) {
        console.log('Speech recognition watchdog: restarting due to inactivity')
        try {
          recognitionRef.current.stop()
          setTimeout(() => {
            if (recognitionRef.current && isSupported) {
              recognitionRef.current.start()
              lastTranscriptTimeRef.current = Date.now()
            }
          }, 1000)
        } catch (error) {
          console.error('Watchdog restart failed:', error)
        }
      }
    }, 10000) // Check every 10 seconds

    return () => clearInterval(watchdogTimer)
  }, [isListening, isSupported])

  return {
    isListening,
    transcript,
    confidence,
    error,
    isSupported,
    startListening,
    stopListening,
    resetTranscript
  }
} 