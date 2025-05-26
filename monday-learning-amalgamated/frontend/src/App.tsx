import React, { useEffect, useState, useCallback } from 'react'
import { Canvas } from '@react-three/fiber'
import { XR, Controllers, Hands } from '@react-three/xr'
import MondayScene from './components/MondayScene'
import VoiceInterface from './components/VoiceInterface'
import SpatialOrchestrator from './components/SpatialOrchestrator'
import LoadingOverlay from './components/LoadingOverlay'
import ErrorBoundary from './components/ErrorBoundary'
import DiagnosticOverlay from './components/DiagnosticOverlay'
import { useMondayStore } from './store/mondayStore'
import { useVoiceSystem } from './hooks/useVoiceSystem'
import { useWebSocketConnection } from './hooks/useWebSocket'
import { usePerformanceMonitor } from './hooks/usePerformanceMonitor'
import { SystemState } from './controllers/VoiceSystemController'

const App: React.FC = () => {
  const [isInitialized, setIsInitialized] = useState(false)
  const [isVRSupported, setIsVRSupported] = useState(false)
  const [lastProcessedTranscript, setLastProcessedTranscript] = useState('')
  const [audioInitialized, setAudioInitialized] = useState(false)
  const [currentModel, setCurrentModel] = useState<string>('')

  const { 
    isConnected, 
    sessionState, 
    initializeSession,
    setPerformanceMetrics,
    addPanel,
    setActivePanel,
    setConversationActive,
    panels
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

  // Initialize the unified voice system
  const {
    systemState,
    isListening,
    isSpeaking,
    isPlaying,
    transcript,
    conversationActive,
    error: voiceError,
    systemStatus,
    handleCommand,
    handleTTSResponse,
    emergencyReset,
    forceStartListening,
    forceStop
  } = useVoiceSystem({
    onError: (error) => {
      console.error('🌟 App: Voice system error:', error)
    }
  })

  // Initialize audio on first user interaction
  const handleUserInteraction = useCallback(async (): Promise<boolean> => {
    if (!audioInitialized) {
      try {
        // Create a temporary audio context to test
        const tempContext = new (window.AudioContext || (window as any).webkitAudioContext)()
        if (tempContext.state === 'suspended') {
          await tempContext.resume()
        }
        await tempContext.close()
        
        setAudioInitialized(true)
        console.log('🌟 App: Audio initialized after user interaction')
        return true
      } catch (error) {
        console.warn('🌟 App: Audio initialization failed:', error)
        return false
      }
    }
    return true
  }, [audioInitialized])

  // Handle interrupting TTS
  const handleInterruptTTS = useCallback(async () => {
    if (systemState === SystemState.PLAYING_TTS) {
      console.log('🌟 App: Interrupting TTS and starting listening')
      await forceStop() // Stop TTS
      await new Promise(resolve => setTimeout(resolve, 500)) // Brief delay
      await forceStartListening() // Start listening
    } else {
      await forceStartListening()
    }
  }, [systemState, forceStop, forceStartListening])

  // Handle WebSocket responses from backend
  useEffect(() => {
    if (!socket) return

    const handleVoiceResponse = async (response: any) => {
      console.log('🌟 App: Received voice response from backend:', response)
      
      // Update current model indicator
      if (response.data?.model) {
        setCurrentModel(response.data.model)
      }
      
      // Add panels if provided
      if (response.data?.panels && Array.isArray(response.data.panels)) {
        console.log('🌟 App: Adding panels to store:', response.data.panels.length, 'panels')
        response.data.panels.forEach((panelData: any, index: number) => {
          console.log(`🌟 App: Adding panel ${index + 1}:`, panelData)
          addPanel(panelData)
        })
        
        const mainPanel = response.data.panels.find((p: any) => p.isActive)
        if (mainPanel) {
          console.log('🌟 App: Setting active panel:', mainPanel.id)
          setActivePanel(mainPanel.id)
        }
      } else {
        console.log('🌟 App: No panels in response or panels is not an array')
      }
      
      // Handle TTS response through the unified system
      if (response.message) {
        console.log('🌟 App: Triggering TTS through voice system controller')
        try {
          await handleTTSResponse(response.message)
          console.log('🌟 App: TTS handling completed successfully')
        } catch (error) {
          console.error('🌟 App: TTS handling failed:', error)
        }
      }
    }

    const handleVoiceError = (errorData: any) => {
      console.error('🌟 App: Backend voice error:', errorData)
    }

    socket.on('voice_response', handleVoiceResponse)
    socket.on('voice_error', handleVoiceError)

    return () => {
      socket.off('voice_response', handleVoiceResponse)
      socket.off('voice_error', handleVoiceError)
    }
  }, [socket, handleTTSResponse, addPanel, setActivePanel])

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

  // Process voice commands when transcript changes (avoids stale closure issues)
  useEffect(() => {
    if (!transcript || transcript === lastProcessedTranscript) {
      return
    }

    console.log('🌟 App: New transcript detected, processing command:', {
      transcript: transcript,
      lastProcessed: lastProcessedTranscript,
      socketConnected,
      conversationActive,
      hasSocket: !!socket,
      socketId: socket?.id || 'no-socket'
    })

    // Check prerequisites with current values (not stale closures)
    if (!socket || !socketConnected || !transcript.trim()) {
      console.log('🌟 App: Command processing blocked - missing prerequisites:', {
        hasSocket: !!socket,
        socketConnected,
        hasCommand: !!transcript.trim(),
        socketReadyState: socket?.connected
      })
      return
    }

    const normalizedCommand = transcript.toLowerCase().trim()
    
    // Check if this should trigger Monday (either explicit trigger or in conversation)
    const isExplicitTrigger = normalizedCommand.includes('hey monday')
    const shouldProcess = isExplicitTrigger || conversationActive
    
    console.log('🌟 App: Command filtering:', {
      command: normalizedCommand.substring(0, 50),
      isExplicitTrigger,
      conversationActive,
      shouldProcess
    })

    if (shouldProcess) {
      // If it's an explicit trigger, set conversation active
      if (isExplicitTrigger) {
        setConversationActive(true)
      }

      console.log('🌟 App: 🎯 Processing voice command:', {
        command: transcript,
        isExplicitTrigger,
        conversationActive,
        commandLength: transcript.length
      })
      
      console.log('🌟 App: 📤 Sending command to backend via WebSocket')
      
      try {
        socket.emit('voice_command', {
          command: transcript,
          timestamp: Date.now(),
          conversationActive: true, // Always send true if we're processing
          isExplicitTrigger: isExplicitTrigger
        })
        
        setLastProcessedTranscript(transcript)
        console.log('🌟 App: ✅ Command sent to backend successfully')
      } catch (error) {
        console.error('🌟 App: ❌ Failed to send command to backend:', error)
      }
    } else {
      console.log('🌟 App: ❌ Command ignored - no trigger and not in conversation:', {
        command: normalizedCommand.substring(0, 50),
        needsTrigger: !isExplicitTrigger && !conversationActive
      })
    }
  }, [transcript, lastProcessedTranscript, socket, socketConnected, conversationActive])

  if (!isInitialized) {
    return <LoadingOverlay message="Initializing Monday..." />
  }

  // Check for critical voice errors that require full-screen handling
  if (voiceError && voiceError.includes('Microphone permission lost')) {
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
          frameloop="always"
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
          onStartListening={handleInterruptTTS}
          onStopListening={forceStop}
          conversationActive={conversationActive}
          onUserInteraction={handleUserInteraction}
          isRestartingVoice={systemState === SystemState.RESETTING}
          isSpeaking={isSpeaking}
          isPlaying={isPlaying}
        />

        {/* Model Indicator */}
        {currentModel && (
          <div style={{
            position: 'fixed',
            top: '1rem',
            right: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'rgba(32, 128, 141, 0.9)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            border: '1px solid var(--true-turquoise)'
          }}>
            🤖 {currentModel.replace('sonar-', '').replace('-', ' ').toUpperCase()}
          </div>
        )}

        {/* Diagnostic Overlay - Only on localhost */}
        {window.location.hostname === 'localhost' && (
          <DiagnosticOverlay
            systemState={systemState}
            systemStatus={systemStatus}
            transcript={transcript}
            conversationActive={conversationActive}
            error={voiceError}
            onEmergencyReset={emergencyReset}
            onForceStart={forceStartListening}
            onForceStop={forceStop}
          />
        )}

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

        {/* Debug: Panel Count */}
        {window.location.hostname === 'localhost' && (
          <div style={{
            position: 'fixed',
            bottom: '1rem',
            right: '8rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'rgba(32, 128, 141, 0.8)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000
          }}>
            Panels: {panels.length}
          </div>
        )}

        {/* Debug: Panel Details */}
        {window.location.hostname === 'localhost' && panels.length > 0 && (
          <div style={{
            position: 'fixed',
            top: '1rem',
            left: '1rem',
            padding: '1rem',
            backgroundColor: 'rgba(9, 23, 23, 0.9)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.75rem',
            zIndex: 1000,
            maxWidth: '300px',
            maxHeight: '200px',
            overflow: 'auto'
          }}>
            <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
              Debug: Panels in Store ({panels.length})
            </div>
            {panels.map((panel, index) => (
              <div key={panel.id} style={{ marginBottom: '0.5rem', padding: '0.25rem', backgroundColor: 'rgba(32, 128, 141, 0.2)' }}>
                <div>#{index + 1}: {panel.title}</div>
                <div>Type: {panel.type}</div>
                <div>ID: {panel.id}</div>
                <div>Active: {panel.isActive ? 'Yes' : 'No'}</div>
                <div>Position: [{panel.position.join(', ')}]</div>
              </div>
            ))}
          </div>
        )}

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

        {/* System Status */}
        <div style={{
          position: 'fixed',
          bottom: '3.5rem',
          right: '1rem',
          padding: '0.5rem 1rem',
          backgroundColor: systemState === SystemState.ERROR ? '#dc3545' : 
                          systemState === SystemState.ACTIVE_LISTENING ? 'var(--true-turquoise)' :
                          systemState === SystemState.PLAYING_TTS ? '#ff6600' :
                          'rgba(9, 23, 23, 0.8)',
          color: 'var(--paper-white)',
          borderRadius: '0.25rem',
          fontSize: '0.875rem',
          zIndex: 1000,
          display: 'flex',
          alignItems: 'center',
          gap: '0.5rem'
        }}>
          {systemState === SystemState.PLAYING_TTS && (
            <div style={{
              width: '8px',
              height: '8px',
              backgroundColor: 'var(--paper-white)',
              borderRadius: '50%',
              animation: 'pulse 1s ease-in-out infinite'
            }} />
          )}
          {systemState === SystemState.ACTIVE_LISTENING ? '🎤 Listening' :
           systemState === SystemState.PLAYING_TTS ? '🔊 Speaking' :
           systemState === SystemState.PROCESSING_COMMAND ? '⚙️ Processing' :
           systemState === SystemState.RESETTING ? '🔄 Resetting' :
           systemState === SystemState.ERROR ? '❌ Error' :
           systemState === SystemState.WAITING_FOR_ACTIVATION ? '⏳ Waiting' :
           '💤 Idle'}
        </div>

        {/* Audio Initialization Prompt */}
        {!audioInitialized && (
          <div 
            onClick={handleUserInteraction}
            style={{
              position: 'fixed',
              bottom: '6rem',
              right: '1rem',
              padding: '1rem',
              backgroundColor: '#f39c12',
              color: 'var(--paper-white)',
              borderRadius: '0.5rem',
              fontSize: '0.875rem',
              zIndex: 1000,
              cursor: 'pointer',
              border: '2px solid #e67e22'
            }}
          >
            🔇 Click to initialize audio for voice responses
          </div>
        )}

        {/* Error Display */}
        {voiceError && !voiceError.includes('Microphone permission lost') && (
          <div style={{
            position: 'fixed',
            bottom: '8rem',
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
          onClick={emergencyReset}
          >
            <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
               Voice System Error
            </div>
            <div style={{ fontSize: '0.75rem', marginBottom: '0.5rem' }}>
              {voiceError}
            </div>
            <div style={{ fontSize: '0.7rem', opacity: 0.8 }}>
              Click to reset the voice system.
            </div>
          </div>
        )}

        {/* Welcome Instructions */}
        {!sessionState.hasCompletedIntro && !conversationActive && systemState === SystemState.WAITING_FOR_ACTIVATION && (
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
            <div style={{
              fontSize: '0.75rem',
              marginTop: '1rem',
              padding: '0.5rem',
              backgroundColor: 'rgba(32, 128, 141, 0.1)',
              borderRadius: '0.25rem'
            }}>
              <strong>System Status:</strong> {systemState} | 
              <strong> Voice:</strong> {isListening ? ' 🎤 Ready' : ' 🔇 Not Active'}
            </div>
          </div>
        )}

        {/* Global Watchdog Alert */}
        {window.location.hostname === 'localhost' && systemStatus.errorCount > 3 && (
          <div 
            onClick={emergencyReset}
            style={{
              position: 'fixed',
              top: '50%',
              right: '1rem',
              transform: 'translateY(-50%)',
              padding: '1rem',
              backgroundColor: '#ff4500',
              color: 'var(--paper-white)',
              borderRadius: '0.5rem',
              fontSize: '0.875rem',
              zIndex: 2001,
              cursor: 'pointer',
              border: '3px solid #ff6500',
              maxWidth: '300px',
              fontWeight: 'bold'
            }}
          >
            🚨 SYSTEM UNSTABLE
            <div style={{ fontSize: '0.7rem', marginTop: '0.5rem', fontWeight: 'normal' }}>
              {systemStatus.errorCount} errors detected. Click for emergency reset.
            </div>
          </div>
        )}
      </div>
    </ErrorBoundary>
  )
}

export default App 