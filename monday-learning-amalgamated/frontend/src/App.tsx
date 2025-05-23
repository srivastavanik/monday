import React, { useEffect, useState } from 'react'
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

const App: React.FC = () => {
  const [isInitialized, setIsInitialized] = useState(false)
  const [isVRSupported, setIsVRSupported] = useState(false)
  const { 
    isConnected, 
    sessionState, 
    initializeSession,
    setPerformanceMetrics 
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
    error: voiceError
  } = useVoiceRecognition()

  // Performance monitoring for Quest optimization
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

  // Auto-start voice recognition when session is ready
  useEffect(() => {
    if (sessionState.isActive && !isListening && !voiceError) {
      startListening()
    }
  }, [sessionState.isActive, isListening, voiceError, startListening])

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
            frameloop="always"
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

        {/* Performance Monitor (Development Only) */}
        {process.env.NODE_ENV === 'development' && (
          <div style={{
            position: 'fixed',
            top: '1rem',
            left: '1rem',
            padding: '0.5rem',
            backgroundColor: 'rgba(9, 23, 23, 0.8)',
            color: 'var(--paper-white)',
            fontSize: '0.75rem',
            borderRadius: '0.25rem',
            fontFamily: 'monospace',
            zIndex: 1000
          }}>
            <div>FPS: {fps}</div>
            <div>Frame: {frameTime.toFixed(2)}ms</div>
            <div>Memory: {memoryUsage.toFixed(1)}MB</div>
            <div>VR: {isVRSupported ? 'Yes' : 'No'}</div>
          </div>
        )}

        {/* Instructions for first-time users */}
        {!sessionState.hasCompletedIntro && (
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
            maxWidth: '400px',
            zIndex: 2000
          }}>
            <h3 style={{ color: 'var(--true-turquoise)', marginBottom: '1rem' }}>
              Welcome to Monday
            </h3>
            <p>Put on your headset and say "Monday, hello" to begin your learning journey.</p>
            <p style={{ fontSize: '0.875rem', opacity: 0.8 }}>
              Make sure microphone access is enabled.
            </p>
          </div>
        )}
      </div>
    </ErrorBoundary>
  )
}

export default App 