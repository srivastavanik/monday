import { useState, useEffect, useCallback, useRef } from 'react'

interface VoiceRecognitionHookProps {
  onListeningStateChange?: (isListening: boolean) => void;
}

interface VoiceRecognitionHook {
  isListening: boolean
  transcript: string
  confidence: number
  error: string | null
  isSupported: boolean
  startListening: () => void
  stopListening: () => void
  resetTranscript: () => void
  clearError: () => void
}

export const useVoiceRecognition = (props?: VoiceRecognitionHookProps): VoiceRecognitionHook => {
  const [isListening, setIsListeningInternal] = useState(false)
  const [transcript, setTranscript] = useState('')
  const [confidence, setConfidence] = useState(0)
  const [error, setError] = useState<string | null>(null)
  const [isSupported, setIsSupported] = useState(false)
  
  const recognitionRef = useRef<any>(null)
  const lastTranscriptTimeRef = useRef<number>(Date.now())
  const isStoppingRef = useRef<boolean>(false)
  const isStartingRef = useRef<boolean>(false)
  const startTimeoutRef = useRef<number | null>(null)
  const globalLockRef = useRef<boolean>(false)

  // Wrapper for setIsListening to call the callback
  const setIsListening = useCallback((value: boolean) => {
    setIsListeningInternal(value);
    if (props?.onListeningStateChange) {
      console.log(`🎤 VoiceRecognition Hook: Calling onListeningStateChange with ${value}`)
      props.onListeningStateChange(value);
    }
  }, [props?.onListeningStateChange]);

  useEffect(() => {
    const SpeechRecognition = (window as any).SpeechRecognition || (window as any).webkitSpeechRecognition
    if (SpeechRecognition) {
      setIsSupported(true)
      recognitionRef.current = new SpeechRecognition()
      console.log('🎤 VoiceRecognition Hook: SpeechRecognition initialized')
    } else {
      setIsSupported(false)
      setError('Speech recognition not supported in this browser')
      console.error('🎤 VoiceRecognition Hook: SpeechRecognition NOT SUPPORTED')
    }
  }, [])

  const startListening = useCallback(() => {
    console.log('🎤 VoiceRecognition Hook: Attempting to call startListening...')
    console.log('🎤 VoiceRecognition Hook: Current state before startListening:', { 
      isListening, // This now refers to the internal state via the getter from useState
      isSupported,
      recognitionRefCurrent: !!recognitionRef.current,
      globalLock: globalLockRef.current,
      isStarting: isStartingRef.current,
      isStopping: isStoppingRef.current,
      startTimeout: startTimeoutRef.current
    })

    if (globalLockRef.current) {
      console.warn('🎤 VoiceRecognition Hook: BLOCKED startListening - Global lock active.')
      return
    }

    if (startTimeoutRef.current) {
      console.log('🎤 VoiceRecognition Hook: Clearing existing startTimeoutRef:', startTimeoutRef.current)
      clearTimeout(startTimeoutRef.current)
      startTimeoutRef.current = null
    }

    if (!recognitionRef.current || !isSupported) {
      console.error('🎤 VoiceRecognition Hook: CANNOT start - Recognition not supported or not initialized.')
      return
    }

    if (isListening) {
      console.warn('🎤 VoiceRecognition Hook: CANNOT start - Already listening (isListening === true).')
      return
    }

    if (isStartingRef.current) {
      console.warn('🎤 VoiceRecognition Hook: CANNOT start - Already starting (isStartingRef.current === true).')
      return
    }

    console.log('🎤 VoiceRecognition Hook: Proceeding with startListening logic...')
    globalLockRef.current = true
    isStartingRef.current = true
    isStoppingRef.current = false 
    console.log('🎤 VoiceRecognition Hook: SET globalLockRef=true, isStartingRef=true, isStoppingRef=false')

    try {
      setError(null)
      setTranscript('') 
      lastTranscriptTimeRef.current = Date.now()
      
      console.log('🎤 VoiceRecognition Hook: Setting isListening to true (PRE-START)')
      setIsListening(true) 
      
      recognitionRef.current.start()
      console.log('🎤 VoiceRecognition Hook: recognition.start() CALLED.')
      
      startTimeoutRef.current = setTimeout(() => {
        if (isStartingRef.current) { 
          console.error('🎤 VoiceRecognition Hook: TIMEOUT - recognition.onstart did not fire in 3s. Resetting state.')
          isStartingRef.current = false
          globalLockRef.current = false
          setIsListening(false) 
          console.log('🎤 VoiceRecognition Hook: RESET isStartingRef=false, globalLockRef=false. SET isListening=false (TIMEOUT)')
          setError('Voice recognition failed to start (timeout).')
        }
        startTimeoutRef.current = null
      }, 3000) as any 
      console.log('🎤 VoiceRecognition Hook: Set startTimeoutRef:', startTimeoutRef.current)

    } catch (err: any) {
      console.error('🎤 VoiceRecognition Hook: ERROR in startListening catch block:', err)
      
      setIsListening(false) 
      isStartingRef.current = false
      globalLockRef.current = false 
      console.log('🎤 VoiceRecognition Hook: RESET isStartingRef=false, globalLockRef=false. SET isListening=false (CATCH ERROR)')
      
      if (err.name === 'InvalidStateError') {
        console.warn('🎤 VoiceRecognition Hook: InvalidStateError - Attempting recovery.')
        setError('Voice recognition state error - auto-recovering')
        
        try {
          if (recognitionRef.current) {
            recognitionRef.current.stop()
            console.log('🎤 VoiceRecognition Hook: Called recognition.stop() for InvalidStateError recovery.')
          }
        } catch (stopErr) {
          console.error('🎤 VoiceRecognition Hook: Error calling stop() during InvalidStateError recovery:', stopErr)
        }
        setTimeout(() => {
          isStartingRef.current = false
          globalLockRef.current = false
          setIsListening(false) 
          console.log('🎤 VoiceRecognition Hook: RESET isStartingRef=false, globalLockRef=false. SET isListening=false (POST InvalidStateError TIMEOUT)')
        }, 500)
      } else {
        setError(`Failed to start voice recognition: ${err.message}`)
      }
    }
  }, [isListening, isSupported, setIsListening]) 

  const stopListening = useCallback(() => {
    console.log('🎤 VoiceRecognition Hook: Attempting to call stopListening...')
    console.log('🎤 VoiceRecognition Hook: Current state before stopListening:', { 
      isListening, 
      isStarting: isStartingRef.current, 
      isStopping: isStoppingRef.current,
      globalLock: globalLockRef.current,
      startTimeout: startTimeoutRef.current 
    })

    if (startTimeoutRef.current) {
      console.log('🎤 VoiceRecognition Hook: Clearing startTimeoutRef in stopListening:', startTimeoutRef.current)
      clearTimeout(startTimeoutRef.current)
      startTimeoutRef.current = null
    }

    if (!recognitionRef.current) {
      console.warn('🎤 VoiceRecognition Hook: CANNOT stop - recognitionRef is null.')
      isStoppingRef.current = false
      isStartingRef.current = false
      globalLockRef.current = false
      if (isListening) {
        setIsListening(false)
        console.log('🎤 VoiceRecognition Hook: SET isListening=false (recognitionRef null but was listening)')
      }
      return
    }

    if (isListening || isStartingRef.current) {
      try {
        console.log('🎤 VoiceRecognition Hook: Proceeding with stopListening logic...')
        isStoppingRef.current = true 
        isStartingRef.current = false 
        console.log('🎤 VoiceRecognition Hook: SET isStoppingRef=true, isStartingRef=false')
        
        recognitionRef.current.stop()
        console.log('🎤 VoiceRecognition Hook: recognition.stop() CALLED.')
        
        if (isListening) { 
            setIsListening(false)
            console.log('🎤 VoiceRecognition Hook: SET isListening=false (IMMEDIATE in stopListening)')
        }
      } catch (error) {
        console.error('🎤 VoiceRecognition Hook: ERROR in stopListening catch block:', error)
        isStoppingRef.current = false
        isStartingRef.current = false 
        globalLockRef.current = false 
        setIsListening(false) 
        console.log('🎤 VoiceRecognition Hook: RESET isStoppingRef=false, isStartingRef=false, globalLockRef=false. SET isListening=false (CATCH ERROR in stop)')
      }
    } else {
      console.log('🎤 VoiceRecognition Hook: Not listening or starting, so just ensuring flags are reset.')
      isStoppingRef.current = false
      isStartingRef.current = false
    }
  }, [isListening, setIsListening]) 

  const resetTranscript = useCallback(() => {
    console.log('🎤 VoiceRecognition Hook: resetTranscript called.')
    setTranscript('')
    setConfidence(0)
  }, [])

  const clearError = useCallback(() => {
    console.log('🎤 VoiceRecognition Hook: clearError called. Current error:', error)
    setError(null)
    
    if (startTimeoutRef.current) {
      clearTimeout(startTimeoutRef.current)
      startTimeoutRef.current = null
      console.log('🎤 VoiceRecognition Hook: Cleared startTimeoutRef in clearError.')
    }
    console.log('🎤 VoiceRecognition Hook: Error cleared.')
  }, [error]) 

  useEffect(() => {
    if (!recognitionRef.current) {
      console.log('🎤 VoiceRecognition Hook: Recognition setup skipped - ref not ready.')
      return
    }
    console.log('🎤 VoiceRecognition Hook: Setting up SpeechRecognition event handlers...')

    const recognition = recognitionRef.current

    recognition.continuous = true
    recognition.interimResults = true
    recognition.lang = 'en-US'

    const onStartHandler = () => {
      console.log('🎤 VoiceRecognition Hook: EVENT onstart FIRED.')
      if (!isListening) { 
          setIsListening(true)
          console.log('🎤 VoiceRecognition Hook: SET isListening=true (IN ONSTART - was false)')
      } else {
          console.log('🎤 VoiceRecognition Hook: isListening already true (IN ONSTART)')
      }
      setError(null) 
      
      isStartingRef.current = false 
      globalLockRef.current = false 
      isStoppingRef.current = false 
      console.log('🎤 VoiceRecognition Hook: RESET isStartingRef=false, globalLockRef=false, isStoppingRef=false (ONSTART)')
      
      if (startTimeoutRef.current) {
        console.log('🎤 VoiceRecognition Hook: Clearing startTimeoutRef in onstart:', startTimeoutRef.current)
        clearTimeout(startTimeoutRef.current)
        startTimeoutRef.current = null
      }
      console.log('🎤 VoiceRecognition Hook: Speech recognition fully started.')
    }

    const onResultHandler = (event: any) => {
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
      if (transcript !== fullTranscript) { 
        setTranscript(fullTranscript)
      }
      
      if (fullTranscript.trim()) {
        lastTranscriptTimeRef.current = Date.now()
      }
    }

    const onErrorHandler = (event: any) => {
      console.error('🎤 VoiceRecognition Hook: EVENT onerror FIRED. Error:', event.error)
      
      if (isStoppingRef.current && (event.error === 'aborted' || event.error === 'no-speech')) { 
        console.log('🎤 VoiceRecognition Hook: Recognition intentionally stopped or aborted while stopping, ignoring error:', event.error)
      } else if (event.error === 'network') {
        setError('Network error in speech recognition. Please check your connection.')
        console.warn('🎤 VoiceRecognition Hook: Network error.')
      } else if (event.error === 'not-allowed' || event.error === 'service-not-allowed') {
        setError('Microphone access denied or speech service not allowed.')
        console.error('🎤 VoiceRecognition Hook: Microphone access denied.')
        setIsSupported(false) 
      } else if (event.error === 'no-speech') {
        setError('No speech detected. Please try speaking louder or closer to the microphone.')
        console.warn('🎤 VoiceRecognition Hook: No speech detected.')
      } else {
        setError(`Speech recognition error: ${event.error}`)
        console.warn('🎤 VoiceRecognition Hook: Other speech error:', event.error)
      }
      
      if (isStartingRef.current) isStartingRef.current = false
      if (globalLockRef.current) globalLockRef.current = false 
      if (event.error === 'not-allowed' || event.error === 'service-not-allowed') {
          if (isListening) setIsListening(false)
          console.log('🎤 VoiceRecognition Hook: SET isListening=false (CRITICAL ERROR like not-allowed)')
      }
      
      console.log('🎤 VoiceRecognition Hook: State after onerror:', {
        isListening, 
        isStarting: isStartingRef.current,
        globalLock: globalLockRef.current,
        isStopping: isStoppingRef.current
      })
       if (startTimeoutRef.current) { 
        clearTimeout(startTimeoutRef.current);
        startTimeoutRef.current = null;
      }
    }

    const onEndHandler = () => {
      console.log('🎤 VoiceRecognition Hook: EVENT onend FIRED.')
      console.log('🎤 VoiceRecognition Hook: State before onend processing:', {
          isListening, 
          isStartingRef: isStartingRef.current,
          isStoppingRef: isStoppingRef.current,
          globalLockRef: globalLockRef.current
      })

      const wasStarting = isStartingRef.current
      const wasStopping = isStoppingRef.current
      
      isStartingRef.current = false
      isStoppingRef.current = false
      globalLockRef.current = false 
      
      if (isListening) { 
        setIsListening(false)
        console.log('🎤 VoiceRecognition Hook: SET isListening=false (IN ONEND - was true)')
      } else {
        console.log('🎤 VoiceRecognition Hook: isListening already false (IN ONEND)')
      }
      
      if (startTimeoutRef.current) {
        console.log('🎤 VoiceRecognition Hook: Clearing startTimeoutRef in onend:', startTimeoutRef.current)
        clearTimeout(startTimeoutRef.current)
        startTimeoutRef.current = null
      }
      
      console.log('🎤 VoiceRecognition Hook: Speech recognition fully ended. Final flags:', {
          isStartingRef: isStartingRef.current,
          isStoppingRef: isStoppingRef.current,
          globalLockRef: globalLockRef.current
      })
      if (wasStarting && !wasStopping) {
          console.warn("🎤 VoiceRecognition Hook: Recognition ended possibly before or shortly after starting properly.")
      }
    }

    recognition.onstart = onStartHandler
    recognition.onresult = onResultHandler
    recognition.onerror = onErrorHandler
    recognition.onend = onEndHandler
    
    console.log('🎤 VoiceRecognition Hook: Event handlers SET UP.')

    return () => {
      console.log('🎤 VoiceRecognition Hook: CLEANUP - Removing event handlers.')
      if (recognitionRef.current) {
        recognitionRef.current.onstart = null
        recognitionRef.current.onresult = null
        recognitionRef.current.onerror = null
        recognitionRef.current.onend = null
        console.log('🎤 VoiceRecognition Hook: Event handlers REMOVED.')
      }
      if (startTimeoutRef.current) {
        clearTimeout(startTimeoutRef.current)
        console.log('🎤 VoiceRecognition Hook: Cleared startTimeoutRef in cleanup.')
      }
    }
  }, [isListening, isSupported, transcript, setIsListening]) 

  useEffect(() => {
    if (!isListening || !recognitionRef.current || !isSupported) return

    console.log('🎤 VoiceRecognition Hook: Watchdog ACTIVATED.')
    const watchdogTimer = setInterval(() => {
      const timeSinceLastTranscript = Date.now() - lastTranscriptTimeRef.current
      
      if (timeSinceLastTranscript > 30000) { 
        console.warn('🎤 VoiceRecognition Hook: WATCHDOG - 30s inactivity. Attempting restart.')
        
        if (globalLockRef.current && !isStoppingRef.current) { // Allow if stopping for watchdog
            console.warn("🎤 VoiceRecognition Hook: WATCHDOG - Restart blocked by global lock (and not already stopping for watchdog).")
            return;
        }

        try {
          console.log('🎤 VoiceRecognition Hook: WATCHDOG - Acquiring lock and setting isStoppingRef=true')
          globalLockRef.current = true;
          isStoppingRef.current = true; 
          recognitionRef.current.stop() 
          console.log('🎤 VoiceRecognition Hook: WATCHDOG - Called stop().')

          setTimeout(() => {
            if (recognitionRef.current && isSupported) {
              // globalLockRef should be released by onEnd. isStoppingRef also.
              // If not, there's a deeper issue, but we attempt to clear it before starting.
              console.log('🎤 VoiceRecognition Hook: WATCHDOG - Releasing lock before calling startListening(). Current lock state:', globalLockRef.current)
              globalLockRef.current = false; 
              isStartingRef.current = false; // Ensure starting is false before calling start
              isStoppingRef.current = false; // Ensure stopping is false before calling start
              
              startListening() 
              lastTranscriptTimeRef.current = Date.now()
              console.log('🎤 VoiceRecognition Hook: WATCHDOG - Restart attempt initiated via startListening().')
            } else {
                 globalLockRef.current = false; 
                 isStoppingRef.current = false; // Also reset if cannot restart
                 console.error("🎤 VoiceRecognition Hook: WATCHDOG - Cannot restart, recognition not supported or ref missing.")
            }
          }, 1000) 
        } catch (error) {
          console.error('🎤 VoiceRecognition Hook: WATCHDOG - Restart failed:', error)
          globalLockRef.current = false 
          isStoppingRef.current = false; 
        }
      }
    }, 10000) 

    return () => {
      clearInterval(watchdogTimer)
      console.log('🎤 VoiceRecognition Hook: Watchdog DEACTIVATED.')
    }
  }, [isListening, isSupported, startListening, setIsListening]) 

  return {
    isListening,
    transcript,
    confidence,
    error,
    isSupported,
    startListening,
    stopListening,
    resetTranscript,
    clearError
  }
} 