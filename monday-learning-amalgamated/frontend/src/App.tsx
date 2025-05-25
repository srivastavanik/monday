import React, { useEffect, useState, useCallback, useRef } from 'react'
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
  const [isRestartingVoice, setIsRestartingVoice] = useState(false)
  
  // This ref will hold the LATEST isListening state from the hook IMMEDIATELY
  const voiceListeningStateRef = useRef(false);

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

  // Callback for the hook to update App.tsx's understanding of isListening
  const handleListeningStateChange = useCallback((listening: boolean) => {
    console.log(`📢 App.tsx: handleListeningStateChange CALLED. New listening state: ${listening}`);
    voiceListeningStateRef.current = listening;
    // Optionally, if you still need a React state for isListening in App.tsx for other UI reasons,
    // you could set it here, but be mindful of potential re-renders.
    // For now, we'll primarily use the ref for immediate checks.
  }, []);

  // Initialize voice recognition with the callback
  const {
    isListening, // This is the state from the hook, might have slight delay
    transcript,
    startListening,
    stopListening,
    error: voiceError,
    resetTranscript,
    clearError
  } = useVoiceRecognition({ onListeningStateChange: handleListeningStateChange });

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
      return true
    }
    return false
  }, [initializeAudio, audioInitialized])

  // Process voice commands and send to backend
  const processVoiceCommand = useCallback(async (command: string) => {
    if (!socket || !socketConnected || !command.trim()) return

    const normalizedCommand = command.toLowerCase().trim()
    
    // Check if this should trigger Monday (either explicit trigger or in conversation)
    const isExplicitTrigger = normalizedCommand.includes('hey monday')
    const shouldProcess = isExplicitTrigger || conversationActive
    
    if (shouldProcess) {
      console.log('📢 App.tsx: 🎯 Processing voice command:', {
        command: command,
        isExplicitTrigger,
        conversationActive,
        commandLength: command.length
      })
      
      if (voiceError) {
        console.log('📢 App.tsx: Clearing voice error before processing command');
        clearError()
      }
      stopTTS() // Stop any ongoing TTS
      
      try {
        await handleUserInteraction() 
        console.log('📢 App.tsx: ✅ Audio ready for response');
      } catch (error) {
        console.warn('📢 App.tsx: ⚠️ Audio initialization failed during command processing:', error)
      }
      
      console.log('📢 App.tsx: 📤 Sending command to backend via WebSocket');
      socket.emit('voice_command', {
        command: command,
        timestamp: Date.now(),
        conversationActive: conversationActive,
        isExplicitTrigger: isExplicitTrigger
      })
      
      resetTranscript()
      setLastProcessedTranscript(command)
      console.log('📢 App.tsx: ✅ Command sent and transcript reset');
    } else {
      console.log('📢 App.tsx: ❌ Voice command ignored (not triggered):', {
        command: normalizedCommand.substring(0, 50),
        hasHeyMonday: isExplicitTrigger,
        conversationActive
      });
    }
  }, [socket, socketConnected, resetTranscript, stopTTS, conversationActive, handleUserInteraction, voiceError, clearError])

  // Function to clear TTS errors
  const clearTTSError = useCallback(() => {
    console.log('📢 App.tsx: Clearing TTS error via stopTTS()');
    stopTTS(); // The stop function in useTextToSpeech clears errors
  }, [stopTTS]);

  // Monitor transcript changes and process commands
  useEffect(() => {
    if (transcript && transcript !== lastProcessedTranscript && transcript.length > 3) { // Reduced length for faster response
      const timer = setTimeout(() => {
        console.log(`📢 App.tsx: Transcript changed to "${transcript}", processing.`);
        processVoiceCommand(transcript)
      }, 700); // Reduced delay
      
      return () => clearTimeout(timer)
    }
  }, [transcript, lastProcessedTranscript, processVoiceCommand])

  // Handle WebSocket responses from backend
  useEffect(() => {
    if (!socket) return

    const handleVoiceResponse = async (response: any) => {
      console.log('📢 App.tsx: Received voice response from backend:', response)
      
      if (response.data?.conversationActive !== undefined) {
        console.log('📢 App.tsx: Setting conversationActive to:', response.data.conversationActive)
        setConversationActive(response.data.conversationActive)
      }
      
      if (voiceError) {
        console.log('📢 App.tsx: Clearing voiceError as Monday is responding.')
        clearError() 
      }
      
      // Stop listening using the ref for the most current state check
      if (voiceListeningStateRef.current) {
        console.log('📢 App.tsx: Voice is listening (checked via ref). Attempting to stop for TTS...')
        stopListening() 
        // The hook's onListeningStateChange will update voiceListeningStateRef.current to false
        // We might not need to aggressively poll here anymore if the callback is reliable.
        // For safety, a very short timeout or a useEffect watching the actual hook's isListening state can be used if needed.
      } else {
        console.log('📢 App.tsx: Voice was not listening (checked via ref) when TTS response received.')
      }
      
      response.data.panels?.forEach((panelData: any) => addPanel(panelData));
      const mainPanel = response.data.panels?.find((p: any) => p.isActive);
      if (mainPanel) setActivePanel(mainPanel.id);
      
      let ttsSucceeded = false;
      console.log('📢 App.tsx: Setting TTS status to speaking')
      setTtsStatus('speaking');

      try {
        const audioReady = await handleUserInteraction();
        console.log('📢 App.tsx: Audio initialized for TTS:', audioReady);
        
        if (audioReady && response.message) {
          console.log('📢 App.tsx: Preparing to speak. Message content:', JSON.stringify(response.message)); // Log the message content
          try {
            await speak(response.message);
            console.log('📢 App.tsx: TTS speak() promise resolved.');
            
            // Check if there was actually a TTS error even though promise resolved
            if (ttsError) {
              console.error('📢 App.tsx: TTS hook reported error after speak():', ttsError);
              ttsSucceeded = false;
            } else {
              ttsSucceeded = true;
            }
          } catch (speakError) {
            console.error('📢 App.tsx: TTS speak() promise rejected:', speakError);
            ttsSucceeded = false;
          }
          
          console.log('📢 App.tsx: Waiting for TTS audio playback to complete (isSpeaking/isPlaying to be false). Current state:', { isSpeaking, isPlaying })
          await new Promise<void>(resolve => {
            const interval = setInterval(() => {
              if (!isSpeaking && !isPlaying) {
                clearInterval(interval);
                console.log('📢 App.tsx: Confirmed TTS audio playback completed.')
                resolve();
              }
            }, 50);
            setTimeout(() => { // Timeout for safety
                clearInterval(interval);
                console.warn('📢 App.tsx: TIMEOUT waiting for TTS playback completion.');
                resolve();
            }, 7000); // 7s timeout, adjust as needed
          });

        } else {
          console.warn('📢 App.tsx: Audio not available or no message for TTS, skipping TTS playback.');
          if (!response.message) console.warn('📢 App.tsx: No message content in response for TTS.');
          ttsSucceeded = false; // No audio available = not successful
        }
      } catch (error) {
        console.error('📢 App.tsx: TTS error during speak() call:', error);
        ttsSucceeded = false;
      } finally {
        console.log('📢 App.tsx: TTS operation finished. TTS Success Flag:', ttsSucceeded);
        const newStatus = ttsSucceeded ? 'ready' : (audioInitialized ? (ttsError ? 'failed' : 'no-audio') : 'no-audio');
        setTtsStatus(newStatus);
        console.log('📢 App.tsx: TTS status set to:', newStatus);
        if (ttsError) console.error('📢 App.tsx: TTS hook reported error:', ttsError);
      }
      
      // Restart logic is now primarily driven by useEffect watching ttsStatus
    }

    socket.on('voice_response', handleVoiceResponse)
    // ... (socket.on('voice_error', handleVoiceError))

    return () => {
      socket.off('voice_response', handleVoiceResponse)
      // ... (socket.off('voice_error', handleVoiceError))
    }
  }, [socket, speak, stopListening, voiceError, clearError, audioInitialized, handleUserInteraction, conversationActive, isSpeaking, isPlaying, ttsError, addPanel, setActivePanel])

  // Single, comprehensive voice restart function - SIMPLIFIED
  const restartVoiceRecognition = useCallback(async (reason: string) => {
    console.log(`📢 App.tsx: 🔄 Voice restart requested. Reason: "${reason}". Current state:`, {
      isListening_hook: isListening,
      isListening_ref: voiceListeningStateRef.current,
      audioInitialized,
      conversationActive,
      ttsStatus,
      isRestartingVoice
    });

    if (!audioInitialized) {
      console.warn(`📢 App.tsx: Voice restart "${reason}" BLOCKED - audio not initialized. Attempting init...`);
      const success = await handleUserInteraction();
      if (!success) {
        console.error(`📢 App.tsx: Voice restart "${reason}" FAILED - audio could not be initialized.`);
        return;
      }
      console.log(`📢 App.tsx: ✅ Audio initialized during voice restart for "${reason}".`);
    }

    // Use the ref for the most immediate check
    if (voiceListeningStateRef.current) {
        console.warn(`📢 App.tsx: Voice restart "${reason}" SKIPPED - already listening (ref: ${voiceListeningStateRef.current}).`);
        return;
    }
    if (isRestartingVoice) {
        console.warn(`📢 App.tsx: Voice restart "${reason}" SKIPPED - already restarting.`);
        return;
    }

    console.log(`📢 App.tsx: 🚀 Starting voice restart for "${reason}"`);
    setIsRestartingVoice(true);

    if (voiceError) {
      console.log(`📢 App.tsx: Clearing voiceError ("${voiceError}") before restart.`);
      clearError();
    }
    
    // Brief delay for stability
    await new Promise(resolve => setTimeout(resolve, 100)); 

    console.log(`📢 App.tsx: 📞 Calling startListening() for "${reason}"`);
    startListening();

    // Enhanced verification that voice recognition actually starts
    let verificationAttempts = 0;
    const maxVerificationAttempts = 20; // 2 seconds max
    
    const verifyStart = async () => {
      while (verificationAttempts < maxVerificationAttempts) {
        await new Promise(resolve => setTimeout(resolve, 100));
        verificationAttempts++;
        
        if (voiceListeningStateRef.current) {
          console.log(`📢 App.tsx: ✅ Voice restart "${reason}" CONFIRMED - listening is active after ${verificationAttempts * 100}ms`);
          if (conversationActive) {
            console.log(`📢 App.tsx: 💬 CONVERSATION READY - User can speak without "Hey Monday"`);
          }
          break;
        }
        
        if (voiceError) {
          console.error(`📢 App.tsx: ❌ Voice restart "${reason}" FAILED - error detected: ${voiceError}`);
          break;
        }
      }
      
      if (verificationAttempts >= maxVerificationAttempts && !voiceListeningStateRef.current) {
        console.warn(`📢 App.tsx: ⚠️ Voice restart "${reason}" TIMEOUT - listening not confirmed after 2s`);
      }
    };
    
    // Start verification process
    verifyStart();

    // Reset restart flag after delay
    setTimeout(() => {
        setIsRestartingVoice(false);
        console.log(`📢 App.tsx: 🏁 Voice restart process "${reason}" completed. Final listening state: ${voiceListeningStateRef.current}`);
    }, 2500);

  }, [audioInitialized, handleUserInteraction, voiceError, clearError, startListening, isRestartingVoice, isListening, conversationActive, ttsStatus]);

  // Auto-start voice recognition when session is ready
  useEffect(() => {
    console.log('📢 App.tsx: Auto-start check. Conditions:', {
      sessionActive: sessionState.isActive,
      isListening: voiceListeningStateRef.current, // Use REF for immediate check
      voiceError: !!voiceError,
      isSpeaking,
      isPlaying,
      isRestartingVoice,
      conversationActive,
      audioInitialized
    })
    if (audioInitialized && 
        sessionState.isActive && 
        !voiceListeningStateRef.current && // Use REF
        !voiceError && 
        !isSpeaking && 
        !isPlaying && 
        !isRestartingVoice && 
        !conversationActive) { 
      console.log('📢 App.tsx: 🚀 Auto-starting voice recognition for new session.');
      restartVoiceRecognition('Auto-start new session');
    } else {
      console.log('📢 App.tsx: Auto-start conditions NOT MET.');
    }
  }, [sessionState.isActive, voiceError, isSpeaking, isPlaying, isRestartingVoice, conversationActive, audioInitialized, restartVoiceRecognition]);

  // Main restart mechanism - triggers when TTS completes and conditions are right.
  useEffect(() => {
    console.log('📢 App.tsx: TTS Completion restart check. Conditions:', {
      ttsStatus,
      conversationActive,
      isListening: voiceListeningStateRef.current, // Use REF
      voiceError: !!voiceError,
      isSpeaking,
      isPlaying,
      isRestartingVoice,
      sessionActive: sessionState.isActive,
      audioInitialized
    })

    if (audioInitialized &&
        ttsStatus === 'ready' && 
        conversationActive && 
        !voiceError && 
        !isSpeaking && 
        !isPlaying && 
        !voiceListeningStateRef.current && // Use REF
        !isRestartingVoice &&
        sessionState.isActive) {
      console.log('📢 App.tsx: 🎯 TTS completion conditions met - calling restartVoiceRecognition.');
      
      // Add a small delay to ensure TTS state has fully settled
      const restartTimer = setTimeout(() => {
        // Double-check conditions are still met
        if (ttsStatus === 'ready' && 
            conversationActive && 
            !voiceListeningStateRef.current && 
            !isRestartingVoice) {
          console.log('📢 App.tsx: 🚀 TTS completion restart confirmed - triggering restart.');
          restartVoiceRecognition('TTS completion');
        } else {
          console.log('📢 App.tsx: ❌ TTS completion restart conditions changed, aborting.');
        }
      }, 200); // Small delay to ensure state stability
      
      return () => clearTimeout(restartTimer);
    } else {
        console.log('📢 App.tsx: TTS completion restart conditions NOT MET.');
    }
  }, [ttsStatus, conversationActive, voiceError, isSpeaking, isPlaying, isRestartingVoice, sessionState.isActive, restartVoiceRecognition, audioInitialized]);

  // Performance monitoring
  const { fps, frameTime, memoryUsage } = usePerformanceMonitor()

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

  if (!isInitialized) {
    return <LoadingOverlay message="Initializing Monday..." />
  }

  if (voiceError && voiceError.includes('not-allowed')) {
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
        <Canvas
          camera={{ 
            position: [0, 1.6, 0],
            fov: 75 
          }}
          gl={{ 
            antialias: true,
            alpha: false,
            powerPreference: 'high-performance' 
          }}
          frameloop="demand"
        >
          <XR referenceSpace="local-floor">
            <ambientLight intensity={0.3} color="#20808D" />
            <directionalLight 
              position={[5, 5, 5]} 
              intensity={0.5} 
              color="#FBFAF4" 
            />
            
            <Controllers />
            <Hands />
            <MondayScene />
            <SpatialOrchestrator />
          </XR>
        </Canvas>

        <VoiceInterface 
          isListening={isListening}
          transcript={transcript}
          onStartListening={startListening}
          onStopListening={stopListening}
          conversationActive={conversationActive}
          onUserInteraction={handleUserInteraction}
          isRestartingVoice={isRestartingVoice}
          isSpeaking={isSpeaking}
          isPlaying={isPlaying}
        />

        {/* Connection Status */}
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

        {/* Conversation Status */}
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

        {/* TTS Status */}
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

        {/* TTS Completion Indicator - only show during restart */}
        {isRestartingVoice && (
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
            ✅ Monday finished - restarting voice recognition...
          </div>
        )}

        {/* Audio Status Indicators */}
        {ttsStatus === 'no-audio' && (
          <div 
            onClick={handleUserInteraction}
            style={{
              position: 'fixed',
              bottom: '6rem',
              right: '1rem',
              padding: '0.5rem 1rem',
              backgroundColor: '#f39c12',
              color: 'var(--paper-white)',
              borderRadius: '0.25rem',
              fontSize: '0.875rem',
              zIndex: 1000,
              cursor: 'pointer'
            }}
          >
            🔇 Audio not initialized - click to enable
          </div>
        )}

        {/* TTS API Error Indicator - ElevenLabs specific */}
        {ttsError && ttsError.includes('ElevenLabs API Error') && (
          <div style={{
            position: 'fixed',
            bottom: '6rem',
            right: '1rem',
            padding: '1rem',
            backgroundColor: '#dc3545',
            color: 'var(--paper-white)',
            borderRadius: '0.5rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            maxWidth: '400px',
            cursor: 'pointer'
          }}
          onClick={clearTTSError}
          >
            <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
              🚫 ElevenLabs TTS Error
            </div>
            <div style={{ fontSize: '0.75rem', marginBottom: '0.5rem' }}>
              {ttsError}
            </div>
            {ttsError.includes('quota') || ttsError.includes('credits') ? (
              <div style={{ fontSize: '0.7rem', marginTop: '0.5rem', opacity: 0.9 }}>
                💡 <strong>Solutions:</strong><br/>
                • Check your ElevenLabs account billing<br/>
                • Add more credits to your account<br/>
                • Use shorter text for TTS requests<br/>
                • Monday will continue working without voice responses
              </div>
            ) : (
              <div style={{ fontSize: '0.7rem', marginTop: '0.5rem', opacity: 0.8 }}>
                Voice responses are disabled until this is resolved.
              </div>
            )}
            <div style={{ fontSize: '0.7rem', marginTop: '0.5rem', opacity: 0.6 }}>
              Click to dismiss.
            </div>
          </div>
        )}

        {/* General TTS Failed Indicator */}
        {ttsStatus === 'failed' && (
          <div 
            onClick={clearTTSError}
            style={{
              position: 'fixed',
              bottom: '6rem',
              right: '1rem',
              padding: '0.5rem 1rem',
              backgroundColor: '#dc3545',
              color: 'var(--paper-white)',
              borderRadius: '0.25rem',
              fontSize: '0.875rem',
              zIndex: 1000,
              cursor: 'pointer'
            }}
          >
            🚫 Voice synthesis failed - click to retry
          </div>
        )}

        {/* Voice Error Indicator */}
        {voiceError && !voiceError.includes('not-allowed') && (
          <div style={{
            position: 'fixed',
            bottom: '8rem',
            right: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: '#f39c12',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            cursor: 'pointer'
          }}
          onClick={clearError}
          >
            ⚠️ Voice issue: {voiceError} (click to dismiss)
          </div>
        )}

        {/* Enhanced Development Monitor with Speech Recognition Debug */}
        {window.location.hostname === 'localhost' && (
          <div 
            onClick={handleUserInteraction}
            style={{
              position: 'fixed',
              top: '1rem',
              left: '1rem',
              padding: '0.5rem',
              backgroundColor: 'rgba(9, 23, 23, 0.95)',
              color: 'var(--paper-white)',
              fontSize: '0.75rem',
              borderRadius: '0.25rem',
              fontFamily: 'monospace',
              zIndex: 1000,
              cursor: 'pointer',
              border: isListening ? '2px solid #00ff00' : '1px solid #20808D',
              maxWidth: '300px'
            }}
          >
            <div style={{ fontWeight: 'bold', color: '#00ff88', marginBottom: '0.5rem' }}>
              SPEECH RECOGNITION DEBUG
            </div>
            <div>FPS: {fps}</div>
            <div>Frame: {frameTime.toFixed(2)}ms</div>
            <div>Memory: {memoryUsage.toFixed(1)}MB</div>
            <div>VR: {isVRSupported ? 'Yes' : 'No'}</div>
            <div style={{ color: isSpeaking || isPlaying ? '#ff6666' : '#66ff66' }}>
              TTS: {isSpeaking ? 'Speaking' : isPlaying ? 'Playing' : 'Ready'}
            </div>
            <div style={{ color: conversationActive ? '#66ff66' : '#ffff66' }}>
              Conv: {conversationActive ? 'Active' : 'Waiting'}
            </div>
            <div style={{ color: audioInitialized ? '#66ff66' : '#ff6666' }}>
              Audio: {audioInitialized ? 'Ready' : 'Not Init'}
            </div>
            <div style={{ color: ttsStatus === 'ready' ? '#66ff66' : '#ff6666' }}>
              Status: {ttsStatus}
            </div>
            <div style={{ 
              color: isListening ? '#00ff00' : '#ff6666',
              fontWeight: 'bold',
              backgroundColor: isListening ? 'rgba(0,255,0,0.1)' : 'rgba(255,0,0,0.1)',
              padding: '0.2rem'
            }}>
              🎤 Listen: {isListening ? 'ON' : 'OFF'}
            </div>
            <div style={{ color: isRestartingVoice ? '#ffff66' : '#66ff66' }}>
              Restart: {isRestartingVoice ? 'YES' : 'NO'}
            </div>
            <div style={{ color: voiceError ? '#ff6666' : '#66ff66' }}>
              Error: {voiceError || 'None'}
            </div>
            {transcript && (
              <div style={{ 
                color: '#66ffff', 
                marginTop: '0.5rem',
                padding: '0.2rem',
                backgroundColor: 'rgba(102,255,255,0.1)',
                borderRadius: '0.2rem',
                fontSize: '0.7rem'
              }}>
                Last: "{transcript.substring(0, 50)}{transcript.length > 50 ? '...' : ''}"
              </div>
            )}
            <div style={{ color: '#ffaa66', marginTop: '0.5rem', fontSize: '0.6rem' }}>
              Click to reset system
            </div>
            <div style={{ marginTop: '0.5rem', display: 'flex', gap: '0.5rem', flexWrap: 'wrap' }}>
              <button
                onClick={(e) => {
                  e.stopPropagation()
                  if (isListening) {
                    stopListening()
                  } else {
                    startListening()
                  }
                }}
                style={{
                  fontSize: '0.6rem',
                  padding: '0.2rem 0.4rem',
                  backgroundColor: isListening ? '#ff4444' : '#44ff44',
                  color: '#000',
                  border: 'none',
                  borderRadius: '0.2rem',
                  cursor: 'pointer'
                }}
              >
                {isListening ? '🛑 Stop' : '🎤 Test'}
              </button>
              <button
                onClick={(e) => {
                  e.stopPropagation()
                  resetTranscript()
                }}
                style={{
                  fontSize: '0.6rem',
                  padding: '0.2rem 0.4rem',
                  backgroundColor: '#4444ff',
                  color: '#fff',
                  border: 'none',
                  borderRadius: '0.2rem',
                  cursor: 'pointer'
                }}
              >
                🧹 Clear
              </button>
            </div>
          </div>
        )}

        {/* Welcome Instructions */}
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
            <p><strong>Say "Hey Monday"</strong> to start your learning journey.</p>
            <p style={{ fontSize: '0.875rem', margin: '1rem 0' }}>
              After that, you can ask questions without saying "Hey Monday" each time:
            </p>
            <ul style={{ fontSize: '0.875rem', textAlign: 'left', color: 'var(--true-turquoise)' }}>
              <li>"Hey Monday, tell me about quantum physics"</li>
              <li>"How does it work?" (no Hey Monday needed)</li>
              <li>"Think about machine learning algorithms"</li>
              <li>"Research the latest developments"</li>
            </ul>
            <p style={{ fontSize: '0.75rem', opacity: 0.8, marginTop: '1rem' }}>
              <strong>Important:</strong> For voice responses, click "Start Listening" first or speak 
              "Hey Monday" to initialize audio.
            </p>
          </div>
        )}
      </div>
    </ErrorBoundary>
  )
}

export default App
 