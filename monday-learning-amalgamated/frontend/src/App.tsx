import React, { useEffect, useState, useCallback } from 'react'
import { Canvas } from '@react-three/fiber'
import { XR, Controllers, Hands } from '@react-three/xr'
import MondayScene from './components/MondayScene'
import VoiceInterface from './components/VoiceInterface'
import SpatialOrchestrator from './components/SpatialOrchestrator'
import LoadingOverlay from './components/LoadingOverlay'
import ErrorBoundary from './components/ErrorBoundary'
import { useMondayStore } from './store/mondayStore'
import { useVoiceRecognition } from './hooks/useVoiceRecognition'
import { useWebSocketConnection } from './hooks/useWebSocket'
import { usePerformanceMonitor } from './hooks/usePerformanceMonitor'
import { useTextToSpeech } from './hooks/useTextToSpeech'

const App: React.FC = () => {
  const [isInitialized, setIsInitialized] = useState(false)
  const [isVRSupported, setIsVRSupported] = useState(false)
  const [lastProcessedTranscript, setLastProcessedTranscript] = useState('')
  const [conversationActive, setConversationActive] = useState(false)
  const [ttsStatus, setTtsStatus] = useState<'ready' | 'speaking' | 'failed' | 'no-audio'>('ready')
  const { 
    isConnected, 
    sessionState, 
    initializeSession,
    setPerformanceMetrics,
    addPanel,
    setActivePanel
  } = useMondayStore()

  // Initialize WebXR and check support
  useEffect(() => {
    const checkVRSupport = async () => {
      if (navigator.xr) {
        try {
          const supported = await navigator.xr.isSessionSupported('immersive-vr')
          setIsVRSupported(supported)
          console.log('VR Support:', supported ? 'Available' : 'Not available')
        } catch (error) {
          console.warn('VR support check failed:', error)
          setIsVRSupported(false)
        }
      } else {
        console.warn('WebXR not available')
        setIsVRSupported(false)
      }
      setIsInitialized(true)
    }

    checkVRSupport()
  }, [])

  // Initialize WebSocket connection
  const { socket, isConnected: socketConnected } = useWebSocketConnection()

  // Initialize voice recognition
  const {
    isListening,
    transcript,
    startListening,
    stopListening,
    error: voiceError,
    resetTranscript
  } = useVoiceRecognition()

  // Initialize text-to-speech
  const { 
    speak, 
    isSpeaking, 
    isPlaying,
    stop: stopTTS,
    error: ttsError,
    initializeAudio
  } = useTextToSpeech({
    stability: 0.6,
    similarityBoost: 0.8,
    speed: 1.1
  })

  // Initialize audio on first user interaction
  const [audioInitialized, setAudioInitialized] = useState(false)
  
  const handleUserInteraction = useCallback(async (): Promise<boolean> => {
    if (!audioInitialized) {
      try {
        const success = await initializeAudio()
        if (success) {
          setAudioInitialized(true)
          console.log('Audio initialized after user interaction')
          return true
        }
      } catch (error) {
        console.warn('Audio initialization failed:', error)
      }
    } else {
      // Already initialized
      return true
    }
    return false
  }, [initializeAudio, audioInitialized])

  // Performance monitoring for Quest optimization
  const { fps, frameTime, memoryUsage } = usePerformanceMonitor()

  // Process voice commands and send to backend
  const processVoiceCommand = useCallback(async (command: string) => {
    if (!socket || !socketConnected || !command.trim()) return

    const normalizedCommand = command.toLowerCase().trim()
    
    // Send command if it contains "Monday" or if we're in an active conversation
    if (normalizedCommand.includes('monday') || conversationActive) {
      console.log('Processing voice command:', command)
      
      // Stop any currently playing TTS
      stopTTS()
      
      // Voice recognition counts as user gesture - try to initialize audio
      try {
        const audioReady = await handleUserInteraction()
        console.log('Audio ready for response:', audioReady)
      } catch (error) {
        console.warn('Audio initialization failed during voice command:', error)
      }
      
      // Send command to backend
      socket.emit('voice_command', {
        command: command,
        timestamp: Date.now()
      })
      
      // Reset transcript after processing
      resetTranscript()
      setLastProcessedTranscript(command)
    }
  }, [socket, socketConnected, resetTranscript, stopTTS, conversationActive, handleUserInteraction])

  // Monitor transcript changes and process commands
  useEffect(() => {
    if (transcript && transcript !== lastProcessedTranscript && transcript.length > 5) {
      // Add a small delay to ensure the transcript is complete
      const timer = setTimeout(() => {
        processVoiceCommand(transcript)
      }, 1000)
      
      return () => clearTimeout(timer)
    }
  }, [transcript, lastProcessedTranscript, processVoiceCommand])

  // Handle WebSocket responses from backend
  useEffect(() => {
    if (!socket) return

    const handleVoiceResponse = async (response: any) => {
      console.log('Received voice response:', response)
      
      if (response.data?.conversationActive !== undefined) {
        setConversationActive(response.data.conversationActive)
      }
      
      if (isListening) {
        stopListening()
        console.log('Stopped listening - Monday is responding')
      }
      
      if (response.data?.panels && response.data.panels.length > 0) {
        console.log('Creating spatial panels:', response.data.panels);
        response.data.panels.forEach((panelData: any) => addPanel(panelData));
        const mainPanel = response.data.panels.find((p: any) => p.isActive);
        if (mainPanel) setActivePanel(mainPanel.id);
      }
      
      let ttsSucceeded = false;
      setTtsStatus('speaking'); // Set status to speaking BEFORE attempting TTS

      try {
        const audioReady = await handleUserInteraction();
        console.log('Audio initialization result:', audioReady);
        
        if (audioReady) {
          await speak(response.message);
          console.log('Monday finished speaking (from speak promise)');
          ttsSucceeded = true;
          
          // Wait for isSpeaking and isPlaying to actually become false
          await new Promise<void>(resolve => {
            const interval = setInterval(() => {
              if (!isSpeaking && !isPlaying) {
                clearInterval(interval);
                resolve();
              }
            }, 50); // Check every 50ms
          });
          console.log('isSpeaking & isPlaying are now false');

        } else {
          console.warn('Audio not available, skipping TTS');
        }
      } catch (error) {
        console.error('TTS error during speak:', error);
        ttsSucceeded = false;
      } finally {
        console.log('TTS operation finished. Success:', ttsSucceeded);
        setTtsStatus(ttsSucceeded ? 'ready' : (audioInitialized ? 'failed' : 'no-audio'));
        console.log('TTS status set to:', ttsSucceeded ? 'ready' : (audioInitialized ? 'failed' : 'no-audio'));
      }
      
      // Log response data for debugging
      if (response.data?.metadata) {
        console.log('Response metadata:', response.data.metadata, 'TTS Actual Success:', ttsSucceeded);
      }
    }

    const handleVoiceError = (error: any) => {
      console.error('Voice command error from backend:', error)
      
      // Resume listening after error
      setTimeout(() => {
        if (!isListening && !voiceError) {
          startListening()
          console.log('Resumed listening after backend error')
        }
      }, 1000)
    }

    socket.on('voice_response', handleVoiceResponse)
    socket.on('voice_error', handleVoiceError)

    return () => {
      socket.off('voice_response', handleVoiceResponse)
      socket.off('voice_error', handleVoiceError)
    }
  }, [socket, speak, isListening, stopListening, startListening, voiceError, isSpeaking, isPlaying, addPanel, setActivePanel])

  // Update performance metrics in store
  useEffect(() => {
    setPerformanceMetrics({
      fps,
      frameTime,
      memoryUsage,
      timestamp: Date.now()
    })
  }, [fps, frameTime, memoryUsage, setPerformanceMetrics])

  // Initialize Monday session
  useEffect(() => {
    if (isInitialized && socketConnected) {
      initializeSession({
        isVRSupported,
        userAgent: navigator.userAgent,
        timestamp: Date.now()
      })
    }
  }, [isInitialized, socketConnected, isVRSupported, initializeSession])

  // Auto-start voice recognition when session is ready
  useEffect(() => {
    if (sessionState.isActive && !isListening && !voiceError && !isSpeaking && !isPlaying) {
      startListening()
    }
  }, [sessionState.isActive, isListening, voiceError, startListening, isSpeaking, isPlaying])

  // Monitor TTS status and reset if stuck
  useEffect(() => {
    const statusCheckTimer = setInterval(() => {
      // Watchdog for TTS actually being stuck in isSpeaking or isPlaying for too long
      if ((isSpeaking || isPlaying) && ttsStatus === 'speaking') {
        // This indicates that the useTextToSpeech hook itself might be stuck
        // Check how long it has been in this state (not implemented here, but could be)
        // For now, if it's speaking/playing but the app thinks it's stuck, it's a deeper issue
        // console.warn('TTS hook reporting speaking/playing, but ttsStatus watchdog trying to reset.');
      } else if (ttsStatus === 'speaking' && !isSpeaking && !isPlaying) {
        // This is the original case: app thinks it's speaking, but hook says it's not.
        console.warn('CRITICAL: TTS status desync. App thinks TTS is speaking, but hook is idle. Resetting to ready.');
        setTtsStatus('ready');
      }

      // Auto-clear 'failed' or 'no-audio' states after a timeout
      if (ttsStatus === 'failed' || ttsStatus === 'no-audio') {
        const clearTimer = setTimeout(() => {
          if (ttsStatus === 'failed' || ttsStatus === 'no-audio') { // Check again before clearing
            console.log('Auto-clearing TTS status:', ttsStatus, 'to ready after timeout');
            setTtsStatus('ready');
          }
        }, 7000); // Clear after 7 seconds if still in this state
        return () => clearTimeout(clearTimer);
      }
    }, 1000); // Check every 1 second

    return () => clearInterval(statusCheckTimer);
  }, [ttsStatus, isSpeaking, isPlaying]);

  // Recovery mechanism: restart voice recognition if stuck
  useEffect(() => {
    const recoveryTimer = setInterval(() => {
      // If we should be listening but aren't, and conversation is active, and TTS is ready (not speaking/failed)
      if (conversationActive && 
          !isListening && 
          !voiceError && 
          ttsStatus === 'ready' && // Ensure TTS is not active or in an error state
          sessionState.isActive
      ) {
        console.warn('Recovery: Voice recognition seems stuck. Restarting.');
        startListening();
      }
    }, 3000); // Check every 3 seconds

    return () => clearInterval(recoveryTimer);
  }, [conversationActive, isListening, voiceError, ttsStatus, sessionState.isActive, startListening]);

  // Handle voice recognition restart after TTS completion
  useEffect(() => {
    // Only restart listening if:
    // 1. TTS is ready (not speaking/failed/no-audio)
    // 2. We're not currently listening
    // 3. There's no voice error
    // 4. We're in an active conversation
    // 5. TTS is not currently active (isSpeaking/isPlaying)
    if (ttsStatus === 'ready' && 
        !isListening && 
        !voiceError && 
        conversationActive && 
        !isSpeaking && 
        !isPlaying) {
      
      console.log('TTS completed, restarting voice recognition');
      
      // Small delay to ensure all TTS processes are fully complete
      const restartTimer = setTimeout(() => {
        if (!isListening && !voiceError && !isSpeaking && !isPlaying) {
          startListening();
          console.log('Voice recognition restarted after TTS completion');
        }
      }, 500);

      return () => clearTimeout(restartTimer);
    }
  }, [ttsStatus, isListening, voiceError, conversationActive, isSpeaking, isPlaying, startListening]);

  // Manual reset function for when system gets stuck
  const manualReset = useCallback(async () => {
    console.log('Manual reset triggered')
    
    // Stop all current activities
    stopTTS()
    stopListening()
    
    // Reset all states
    setTtsStatus('ready')
    setConversationActive(false)
    resetTranscript()
    
    // Clear any pending timers by forcing a state reset
    setTimeout(async () => {
      // Try to reinitialize audio
      try {
        await handleUserInteraction()
      } catch (error) {
        console.warn('Audio reinitialization failed during reset:', error)
      }
      
      // Restart voice recognition
      setTimeout(() => {
        if (!voiceError) {
          startListening()
          console.log('Voice recognition restarted after manual reset')
        }
      }, 1000)
    }, 500)
  }, [stopTTS, stopListening, resetTranscript, handleUserInteraction, startListening, voiceError])

  if (!isInitialized) {
    return <LoadingOverlay message="Initializing Monday..." />
  }

  if (voiceError) {
    return (
      <div style={{ 
        padding: '2rem', 
        textAlign: 'center', 
        color: 'var(--paper-white)',
        backgroundColor: 'var(--offblack)',
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center'
      }}>
        <h2>Microphone Access Required</h2>
        <p>Monday needs microphone access to function properly.</p>
        <p>Please refresh the page and allow microphone access when prompted.</p>
        <button 
          className="btn" 
          onClick={() => window.location.reload()}
          style={{ margin: '1rem auto', display: 'block' }}
        >
          Refresh Page
        </button>
      </div>
    )
  }

  return (
    <ErrorBoundary>
      <div style={{ 
        width: '100vw', 
        height: '100vh', 
        backgroundColor: 'var(--offblack)' 
      }}>
        {/* WebXR Canvas */}
        <Canvas
          camera={{ 
            position: [0, 1.6, 0],  // Average eye height
            fov: 75 
          }}
          gl={{ 
            antialias: true,
            alpha: false,
            powerPreference: 'high-performance' 
          }}
          frameloop="demand" // Optimize for Quest
        >
          <XR 
            referenceSpace="local-floor"
          >
            {/* Lighting setup for VR */}
            <ambientLight intensity={0.3} color="#20808D" />
            <directionalLight 
              position={[5, 5, 5]} 
              intensity={0.5} 
              color="#FBFAF4" 
            />
            
            {/* VR Controllers and Hand Tracking */}
            <Controllers />
            <Hands />
            
            {/* Main Monday Scene */}
            <MondayScene />
            
            {/* Spatial Information Manager */}
            <SpatialOrchestrator />
          </XR>
        </Canvas>

        {/* Voice Interface Overlay */}
        <VoiceInterface 
          isListening={isListening}
          transcript={transcript}
          onStartListening={startListening}
          onStopListening={stopListening}
          conversationActive={conversationActive}
          onUserInteraction={handleUserInteraction}
        />

        {/* Connection Status Indicator */}
        <div style={{
          position: 'fixed',
          bottom: '1rem',
          right: '1rem',
          padding: '0.5rem 1rem',
          backgroundColor: socketConnected ? 'var(--true-turquoise)' : '#dc3545',
          color: 'var(--paper-white)',
          borderRadius: '0.25rem',
          fontSize: '0.875rem',
          zIndex: 1000
        }}>
          {socketConnected ? 'Connected' : 'Disconnected'}
        </div>

        {/* Conversation Status Indicator */}
        {conversationActive && (
          <div style={{
            position: 'fixed',
            bottom: '1rem',
            left: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'var(--true-turquoise)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            display: 'flex',
            alignItems: 'center',
            gap: '0.5rem'
          }}>
            <div style={{
              width: '8px',
              height: '8px',
              backgroundColor: 'var(--paper-white)',
              borderRadius: '50%',
              animation: 'pulse 2s ease-in-out infinite'
            }} />
            Conversation Active
          </div>
        )}

        {/* TTS Status Indicator */}
        {(isSpeaking || isPlaying) && (
          <div style={{
            position: 'fixed',
            bottom: '3.5rem',
            right: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'var(--true-turquoise)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            display: 'flex',
            alignItems: 'center',
            gap: '0.5rem'
          }}>
            <div style={{
              width: '8px',
              height: '8px',
              backgroundColor: 'var(--paper-white)',
              borderRadius: '50%',
              animation: 'pulse 1s ease-in-out infinite'
            }} />
            Monday Speaking
          </div>
        )}

        {/* TTS Status Indicators */}
        {ttsStatus === 'no-audio' && (
          <div 
            onClick={handleUserInteraction}
            style={{
              position: 'fixed',
              bottom: '3.5rem',
              right: '1rem',
              padding: '0.5rem 1rem',
              backgroundColor: '#f39c12',
              color: 'var(--paper-white)',
              borderRadius: '0.25rem',
              fontSize: '0.875rem',
              zIndex: 1000,
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem',
              cursor: 'pointer'
            }}
          >
            🔇 Audio not initialized - click to enable
          </div>
        )}

        {ttsStatus === 'failed' && (
          <div 
            onClick={handleUserInteraction}
            style={{
              position: 'fixed',
              bottom: '3.5rem',
              right: '1rem',
              padding: '0.5rem 1rem',
              backgroundColor: '#dc3545',
              color: 'var(--paper-white)',
              borderRadius: '0.25rem',
              fontSize: '0.875rem',
              zIndex: 1000,
              display: 'flex',
              alignItems: 'center',
              gap: '0.5rem',
              cursor: 'pointer'
            }}
          >
            🚫 Voice synthesis failed - click to retry
          </div>
        )}

        {/* TTS Error Indicator */}
        {ttsError && (
          <div style={{
            position: 'fixed',
            bottom: '3.5rem',
            right: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: '#dc3545',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000
          }}>
            TTS Error: {ttsError}
          </div>
        )}

        {/* Manual Reset Button - visible when system might be stuck */}
        {(ttsStatus === 'speaking' && !isSpeaking && !isPlaying) || 
         (conversationActive && !isListening && !isSpeaking && !isPlaying) && (
          <div style={{
            position: 'fixed',
            bottom: '8rem',
            left: '50%',
            transform: 'translateX(-50%)',
            zIndex: 1001
          }}>
            <button
              onClick={manualReset}
              style={{
                padding: '0.75rem 1.5rem',
                backgroundColor: '#dc3545',
                color: 'var(--paper-white)',
                border: 'none',
                borderRadius: '0.5rem',
                fontSize: '0.875rem',
                cursor: 'pointer',
                fontWeight: 'bold'
              }}
            >
              🔄 Reset System
            </button>
          </div>
        )}

        {/* Performance Monitor (Development Only) */}
        {window.location.hostname === 'localhost' && (
          <div 
            onClick={manualReset}
            style={{
              position: 'fixed',
              top: '1rem',
              left: '1rem',
              padding: '0.5rem',
              backgroundColor: 'rgba(9, 23, 23, 0.8)',
              color: 'var(--paper-white)',
              fontSize: '0.75rem',
              borderRadius: '0.25rem',
              fontFamily: 'monospace',
              zIndex: 1000,
              cursor: 'pointer',
              border: ttsStatus === 'speaking' && !isSpeaking && !isPlaying ? '2px solid #dc3545' : '1px solid #20808D'
            }}
            title="Click to reset system if stuck"
          >
            <div>FPS: {fps}</div>
            <div>Frame: {frameTime.toFixed(2)}ms</div>
            <div>Memory: {memoryUsage.toFixed(1)}MB</div>
            <div>VR: {isVRSupported ? 'Yes' : 'No'}</div>
            <div>TTS: {isSpeaking ? 'Speaking' : isPlaying ? 'Playing' : 'Ready'}</div>
            <div>Conv: {conversationActive ? 'Active' : 'Waiting'}</div>
            <div>Audio: {audioInitialized ? 'Ready' : 'Not Init'}</div>
            <div>Status: {ttsStatus}</div>
            <div>Listen: {isListening ? 'ON' : 'OFF'}</div>
            {(ttsStatus === 'speaking' && !isSpeaking && !isPlaying) && (
              <div style={{ color: '#dc3545', fontWeight: 'bold' }}>⚠️ STUCK</div>
            )}
          </div>
        )}

        {/* Instructions for first-time users */}
        {!sessionState.hasCompletedIntro && !conversationActive && (
          <div style={{
            position: 'fixed',
            top: '50%',
            left: '50%',
            transform: 'translate(-50%, -50%)',
            padding: '2rem',
            backgroundColor: 'var(--paper-white)',
            color: 'var(--offblack)',
            borderRadius: '0.5rem',
            border: '2px solid var(--true-turquoise)',
            textAlign: 'center',
            maxWidth: '500px',
            zIndex: 2000
          }}>
            <h3 style={{ color: 'var(--true-turquoise)', marginBottom: '1rem' }}>
              Welcome to Monday
            </h3>
            <p><strong>Say "Monday, hello"</strong> to start your learning journey.</p>
            <p style={{ fontSize: '0.875rem', margin: '1rem 0' }}>
              After that, you can ask questions without saying "Monday" each time:
            </p>
            <ul style={{ fontSize: '0.875rem', textAlign: 'left', color: 'var(--true-turquoise)' }}>
              <li>"Monday, tell me about quantum physics"</li>
              <li>"How does it work?" (no Monday needed)</li>
              <li>"Think about machine learning algorithms"</li>
              <li>"Research the latest developments"</li>
            </ul>
            <p style={{ fontSize: '0.75rem', opacity: 0.8, marginTop: '1rem' }}>
              <strong>Important:</strong> For voice responses, click "Start Listening" first or speak 
              "Monday, hello" to initialize audio.
            </p>
            {ttsError && (
              <p style={{ fontSize: '0.75rem', color: '#dc3545', marginTop: '0.5rem' }}>
                Audio issue: {ttsError}
              </p>
            )}
            {ttsStatus === 'no-audio' && (
              <p style={{ fontSize: '0.75rem', color: '#f39c12', marginTop: '0.5rem' }}>
                ⚠️ Audio not initialized - Monday can hear you but can't speak back yet
              </p>
            )}
          </div>
        )}
      </div>
    </ErrorBoundary>
  )
}

export default App 