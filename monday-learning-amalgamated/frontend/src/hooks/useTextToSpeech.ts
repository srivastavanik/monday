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
  const [isPlaying, setIsPlayingInternal] = useState(false)
  const [isSpeaking, setIsSpeakingInternal] = useState(false)
  const [error, setError] = useState<string | null>(null)
  
  const wsRef = useRef<WebSocket | null>(null)
  const audioContextRef = useRef<AudioContext | null>(null)
  const audioQueueRef = useRef<ArrayBuffer[]>([])
  const isPlayingAudioRef = useRef(false)
  const audioInitializedRef = useRef(false)
  const ttsGenerationCompleteRef = useRef(false)
  const allAudioBuffersRef = useRef<ArrayBuffer[]>([])
  const currentPlaybackPromiseRef = useRef<Promise<void> | null>(null)

  // Default configuration
  const {
    voiceId = (import.meta as any).env?.VITE_ELEVENLABS_VOICE_ID || 'XrExE9yKIg1WjnnlVkGX',
    apiKey = (import.meta as any).env?.VITE_ELEVENLABS_API_KEY || '',
    outputFormat = 'mp3_44100_128',
    stability = 0.5,
    similarityBoost = 0.8,
    speed = 1.0
  } = config

  const setIsPlaying = useCallback((value: boolean) => {
    console.log(`TTS_HOOK: setIsPlaying called with ${value}`);
    setIsPlayingInternal(value);
  }, []);

  const setIsSpeaking = useCallback((value: boolean) => {
    console.log(`TTS_HOOK: setIsSpeaking (API activity) called with ${value}`);
    setIsSpeakingInternal(value);
  }, []);

  const initializeAudio = useCallback(async (): Promise<boolean> => {
    console.log("TTS_HOOK: Attempting to initialize audio...");
    if (audioInitializedRef.current && audioContextRef.current?.state === 'running') {
      console.log("TTS_HOOK: Audio already initialized and running.");
      return true
    }

    try {
      if (!audioContextRef.current) {
        audioContextRef.current = new (window.AudioContext || (window as any).webkitAudioContext)()
        console.log("TTS_HOOK: New AudioContext created.");
      }
      
      if (audioContextRef.current.state === 'suspended') {
        console.log("TTS_HOOK: AudioContext is suspended, attempting to resume...");
        await audioContextRef.current.resume()
        console.log("TTS_HOOK: AudioContext resume attempt finished. New state:", audioContextRef.current.state);
      }
      
      if (audioContextRef.current.state === 'running') {
        audioInitializedRef.current = true
        console.log('TTS_HOOK: AudioContext initialized and running successfully.')
        setError(null)
        return true
      } else {
        console.error(`TTS_HOOK: AudioContext failed to run. Current state: ${audioContextRef.current.state}`)
        throw new Error(`AudioContext state: ${audioContextRef.current.state}`)
      }
    } catch (err) {
      console.error('TTS_HOOK: Failed to initialize audio context:', err)
      setError('Audio initialization failed - please interact with the page first')
      audioInitializedRef.current = false
      return false
    }
  }, [])

  const playAudioBuffer = useCallback(async (audioBuffer: ArrayBuffer): Promise<void> => {
    if (!audioContextRef.current || audioContextRef.current.state !== 'running') {
      console.error('TTS_HOOK: playAudioBuffer - AudioContext not ready or not running.');
      throw new Error('AudioContext not ready')
    }
    console.log('TTS_HOOK: playAudioBuffer - Decoding audio data...');
    try {
      const audioData = await audioContextRef.current.decodeAudioData(audioBuffer.slice(0))
      const source = audioContextRef.current.createBufferSource()
      source.buffer = audioData
      source.connect(audioContextRef.current.destination)
      
      console.log('TTS_HOOK: playAudioBuffer - Starting playback of one buffer.');
      return new Promise<void>((resolve, reject) => {
        source.onended = () => {
          console.log('TTS_HOOK: playAudioBuffer - Playback of one buffer ended.');
          resolve();
        }
        source.start()
      }).catch(err => {
        console.error('TTS_HOOK: playAudioBuffer - Error during source.start() or onended: ', err);
        throw err;
      });
    } catch (err) {
      console.error('TTS_HOOK: playAudioBuffer - Failed to decode or play audio:', err)
      throw err
    }
  }, [])

  const processAudioQueue = useCallback(async () => {
    if (!ttsGenerationCompleteRef.current) {
        console.log("TTS_HOOK: processAudioQueue - Generation not complete, deferring playback.");
        return;
    }
    if (isPlayingAudioRef.current) {
        console.log("TTS_HOOK: processAudioQueue - Already playing audio, deferring.");
        return;
    }
    if (allAudioBuffersRef.current.length === 0) {
        console.log("TTS_HOOK: processAudioQueue - No audio buffers to play.");
        setIsPlaying(false);
        setIsSpeaking(false);
        ttsGenerationCompleteRef.current = false;
        return;
    }

    console.log(`TTS_HOOK: processAudioQueue - Starting full playback of ${allAudioBuffersRef.current.length} collected buffers.`);
    isPlayingAudioRef.current = true
    setIsPlaying(true)

    try {
      for (let i = 0; i < allAudioBuffersRef.current.length; i++) {
        console.log(`TTS_HOOK: processAudioQueue - Playing buffer ${i + 1} of ${allAudioBuffersRef.current.length}`);
        await playAudioBuffer(allAudioBuffersRef.current[i])
      }
      console.log('TTS_HOOK: processAudioQueue - All audio buffers playback completed successfully.')
    } catch (err) {
      console.error('TTS_HOOK: processAudioQueue - Error during playback of audio buffers:', err)
      setError('Audio playback failed during queue processing.')
    } finally {
      console.log('TTS_HOOK: processAudioQueue - Playback attempt finished. Resetting playback flags.');
      isPlayingAudioRef.current = false
      setIsPlaying(false)
      
      setIsSpeaking(false)
      
      console.log('TTS_HOOK: ✅ TTS COMPLETELY FINISHED (generation and playback).')
      
      ttsGenerationCompleteRef.current = false
      allAudioBuffersRef.current = [] 
      audioQueueRef.current = []
    }
  }, [playAudioBuffer, setIsPlaying, setIsSpeaking])

  const base64ToArrayBuffer = useCallback((base64: string): ArrayBuffer => {
    const binaryString = window.atob(base64)
    const len = binaryString.length
    const bytes = new Uint8Array(len)
    for (let i = 0; i < len; i++) {
      bytes[i] = binaryString.charCodeAt(i)
    }
    return bytes.buffer
  }, [])

  const speak = useCallback(async (text: string): Promise<void> => {
    console.log(`TTS_HOOK: speak() called with text: "${text ? text.substring(0,50) + '...' : 'EMPTY TEXT'}"`);
    if (!apiKey) {
      console.error('TTS_HOOK: speak - API key not configured.')
      setError('ElevenLabs API key not configured')
      setIsSpeaking(false)
      return
    }

    if (!text || !text.trim()) {
      console.warn('TTS_HOOK: speak - Received empty or whitespace-only text. Nothing to speak.')
      setIsSpeaking(false)
      return
    }
    
    if (wsRef.current) {
        console.log("TTS_HOOK: speak - Closing existing WebSocket connection before new speak call.");
        wsRef.current.onclose = null;
        wsRef.current.close();
        wsRef.current = null;
    }
    if (isPlayingAudioRef.current) {
        console.log("TTS_HOOK: speak - Stopping ongoing audio playback due to new speak call.");
    }
    audioQueueRef.current = [];
    allAudioBuffersRef.current = [];
    ttsGenerationCompleteRef.current = false;
    setIsPlaying(false);

    console.log("TTS_HOOK: speak - Setting isSpeaking (API) to true.");
    setError(null)
    setIsSpeaking(true)

    try {
      const audioReady = await initializeAudio()
      if (!audioReady) {
        console.warn('TTS_HOOK: speak - Audio not ready after init attempt. TTS will generate but not play.')
      }

      const wsUrl = `wss://api.elevenlabs.io/v1/text-to-speech/${voiceId}/stream-input?output_format=${outputFormat}&auto_mode=true`
      console.log("TTS_HOOK: speak - Connecting to WebSocket:", wsUrl);
      const ws = new WebSocket(wsUrl)
      wsRef.current = ws

      currentPlaybackPromiseRef.current = new Promise<void>((resolve, reject) => {
        let hasFinishedProcessingMessage = false;

        ws.onopen = () => {
          console.log('TTS_HOOK: ElevenLabs WebSocket ONOPEN.')
          ws.send(JSON.stringify({
            text: ' ',
            voice_settings: { stability, similarity_boost: similarityBoost, speed },
            xi_api_key: apiKey
          }))
          ws.send(JSON.stringify({ text: text + ' ', try_trigger_generation: true }))
          ws.send(JSON.stringify({ text: '' }))
          console.log('TTS_HOOK: WebSocket - Initial messages sent (BOS, text, EOS).');
        }

        ws.onmessage = (event) => {
          try {
            const data = JSON.parse(event.data as string)
            
            if (data.audio) {
              const audioBuffer = base64ToArrayBuffer(data.audio)
              allAudioBuffersRef.current.push(audioBuffer)
              console.log(`TTS_HOOK: WebSocket - Audio chunk ${allAudioBuffersRef.current.length} collected (${audioBuffer.byteLength} bytes).`);
            }
            
            if (data.isFinal === true) {
              console.log(`TTS_HOOK: WebSocket - 'isFinal: true' received. Total chunks: ${allAudioBuffersRef.current.length}.`);
              hasFinishedProcessingMessage = true;
              ttsGenerationCompleteRef.current = true
              
              if (!audioReady || allAudioBuffersRef.current.length === 0) {
                console.warn('TTS_HOOK: WebSocket - No audio to play (audio not ready or no chunks received). Finishing TTS early.');
                setIsSpeaking(false);
                setIsPlaying(false); 
                ttsGenerationCompleteRef.current = false;
                if (wsRef.current) wsRef.current = null;
                resolve();
              } else {
                console.log('TTS_HOOK: WebSocket - Triggering processAudioQueue.');
                processAudioQueue().then(resolve).catch(reject);
              }
            } else if (data.isFinal === false) {
            } else if (data.alignment || data.normalizedAlignment) {
            } else {
            }
          } catch (err) {
            console.error('TTS_HOOK: WebSocket - Error processing ONMESSAGE:', err, "Raw data:", event.data)
          }
        }

        ws.onerror = (event) => {
          console.error('TTS_HOOK: ElevenLabs WebSocket ONERROR:', event)
          setError('TTS WebSocket connection error.')
          setIsSpeaking(false)
          setIsPlaying(false)
          ttsGenerationCompleteRef.current = false
          allAudioBuffersRef.current = []
          if (wsRef.current) wsRef.current = null;
          reject(new Error('TTS WebSocket error'))
        }

        ws.onclose = (event) => {
          console.log(`TTS_HOOK: ElevenLabs WebSocket ONCLOSE. Code: ${event.code}, Reason: "${event.reason}", WasClean: ${event.wasClean}, hasFinishedMsg: ${hasFinishedProcessingMessage}`)
          wsRef.current = null
          
          // Check for API errors (quota, billing, etc.)
          if (event.code === 1008 || event.reason.includes('quota') || event.reason.includes('credits') || event.reason.includes('billing')) {
            console.error('TTS_HOOK: ElevenLabs API quota/billing error:', event.reason)
            setError(`ElevenLabs API Error: ${event.reason}`)
            setIsSpeaking(false);
            setIsPlaying(false);
            ttsGenerationCompleteRef.current = false;
            allAudioBuffersRef.current = [];
            reject(new Error(`ElevenLabs API Error: ${event.reason}`));
            return;
          }
          
          // Check for other API errors (codes 1000-1999 are typically server/API issues)
          if (event.code >= 1000 && event.code < 2000 && event.code !== 1000 && event.code !== 1001) {
            console.error('TTS_HOOK: ElevenLabs API error:', event.code, event.reason)
            setError(`ElevenLabs API Error (${event.code}): ${event.reason || 'Unknown API error'}`)
            setIsSpeaking(false);
            setIsPlaying(false);
            ttsGenerationCompleteRef.current = false;
            allAudioBuffersRef.current = [];
            reject(new Error(`ElevenLabs API Error (${event.code}): ${event.reason || 'Unknown API error'}`));
            return;
          }
          
          if (!hasFinishedProcessingMessage) {
            console.warn("TTS_HOOK: WebSocket closed before 'isFinal: true' message was fully processed.");
            if (audioReady && allAudioBuffersRef.current.length > 0) {
              console.log("TTS_HOOK: Attempting to play buffered audio despite early WS close.");
              ttsGenerationCompleteRef.current = true;
              processAudioQueue().then(resolve).catch(reject);
            } else {
              console.warn("TTS_HOOK: No audio to play after early close. Resolving with no audio.");
              setIsSpeaking(false);
              setIsPlaying(false);
              ttsGenerationCompleteRef.current = false;
              allAudioBuffersRef.current = [];
              resolve();
            }
          } else {
            console.log("TTS_HOOK: WebSocket closed, 'isFinal' block was processed. Assuming completion handled.");
          }
        }
      });
      return currentPlaybackPromiseRef.current;

    } catch (err) {
      console.error('TTS_HOOK: speak - Outer catch block error:', err)
      setError('Text-to-speech failed before connection could be established.')
      setIsSpeaking(false)
      setIsPlaying(false)
      throw err
    }
  }, [apiKey, voiceId, outputFormat, stability, similarityBoost, speed, initializeAudio, base64ToArrayBuffer, processAudioQueue, setIsPlaying, setIsSpeaking])

  const stop = useCallback(() => {
    console.log("TTS_HOOK: stop() called.");
    if (wsRef.current) {
      console.log("TTS_HOOK: Closing WebSocket from stop().");
      wsRef.current.onclose = null;
      wsRef.current.close()
      wsRef.current = null
    }
    
    audioQueueRef.current = []
    allAudioBuffersRef.current = []
    isPlayingAudioRef.current = false
    ttsGenerationCompleteRef.current = false
    
    setIsPlaying(false)
    setIsSpeaking(false)
    setError(null)
    console.log("TTS_HOOK: Playback and API flags reset by stop().");
  }, [setIsPlaying, setIsSpeaking])

  return {
    isPlaying,
    isSpeaking,
    speak,
    stop,
    error,
    initializeAudio
  }
} 