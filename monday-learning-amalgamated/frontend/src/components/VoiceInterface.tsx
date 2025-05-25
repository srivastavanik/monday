import React from 'react'

interface VoiceInterfaceProps {
  isListening: boolean
  transcript: string
  onStartListening: () => void
  onStopListening: () => void
  conversationActive?: boolean
  onUserInteraction?: () => void
  isRestartingVoice?: boolean
  isSpeaking?: boolean
  isPlaying?: boolean
}

const VoiceInterface = ({ 
  isListening, 
  transcript, 
  onStartListening, 
  onStopListening,
  conversationActive = false,
  onUserInteraction,
  isRestartingVoice = false,
  isSpeaking = false,
  isPlaying = false
}: VoiceInterfaceProps) => {

  // Smart button logic - prevent conflicts with the atomic state machine
  const isSystemBusy = isRestartingVoice || isSpeaking || isPlaying
  const canManuallyStart = !isListening && !isSystemBusy
  const shouldShowStartButton = canManuallyStart && !conversationActive // Only show when not in conversation
  
  const handleStartListening = () => {
    // Double-check conditions before attempting start
    if (!canManuallyStart) {
      console.log('🎯 VoiceInterface: Manual start blocked - system busy')
      return
    }
    
    onUserInteraction?.()
    onStartListening()
  }

  const handleStopListening = () => {
    onUserInteraction?.()
    onStopListening()
  }

  // Enhanced status message logic
  const getStatusMessage = () => {
    if (isSpeaking) {
      return 'Monday is generating voice response...'
    }
    if (isPlaying) {
      return 'Monday is speaking...'
    }
    if (isRestartingVoice) {
      return 'Restarting voice recognition...'
    }
    if (isListening) {
      if (conversationActive) {
        return 'Listening... (no "Hey Monday" needed)'
      } else {
        return 'Listening for "Hey Monday..."'
      }
    }
    if (conversationActive) {
      return 'In conversation - voice ready to restart'
    }
    return 'Voice system ready'
  }

  const getStatusColor = () => {
    if (isSpeaking || isPlaying) return 'var(--true-turquoise)'
    if (isRestartingVoice) return '#ffaa00'
    if (isListening) return 'var(--true-turquoise)'
    if (conversationActive) return '#00cc66'
    return 'rgba(9, 23, 23, 0.8)'
  }

  return (
    <div style={{
      position: 'fixed',
      bottom: '2rem',
      left: '50%',
      transform: 'translateX(-50%)',
      padding: '1rem 2rem',
      backgroundColor: getStatusColor(),
      color: 'var(--paper-white)',
      borderRadius: '2rem',
      border: `2px solid ${isListening ? 'var(--true-turquoise)' : '#20808D'}`,
      minWidth: '300px',
      textAlign: 'center',
      zIndex: 1000,
      backdropFilter: 'blur(10px)',
      transition: 'all 0.3s ease'
    }}>
      {/* Status Message */}
      <div style={{
        fontSize: '0.875rem',
        marginBottom: '0.5rem',
        opacity: 0.9,
        fontWeight: '500'
      }}>
        {getStatusMessage()}
      </div>
      
      {/* Current Transcript */}
      {transcript && (
        <div style={{
          fontSize: '1rem',
          fontWeight: 500,
          padding: '0.5rem 0',
          backgroundColor: 'rgba(255,255,255,0.1)',
          borderRadius: '0.5rem',
          margin: '0.5rem 0'
        }}>
          "{transcript}"
        </div>
      )}
      
      {/* Conversation Indicator */}
      {conversationActive && (
        <div style={{
          fontSize: '0.75rem',
          color: '#00ff88',
          marginBottom: '0.5rem',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          gap: '0.5rem'
        }}>
          <div style={{
            width: '6px',
            height: '6px',
            backgroundColor: '#00ff88',
            borderRadius: '50%',
            animation: 'pulse 2s ease-in-out infinite'
          }} />
          💬 Continue speaking without "Hey Monday"
        </div>
      )}
      
      {/* Controls */}
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        gap: '1rem',
        marginTop: '0.5rem'
      }}>
        {shouldShowStartButton ? (
          <button
            className="btn"
            onClick={handleStartListening}
            style={{
              fontSize: '0.875rem',
              padding: '0.5rem 1rem',
              backgroundColor: 'var(--paper-white)',
              color: 'var(--true-turquoise)',
              border: '2px solid var(--true-turquoise)',
              borderRadius: '0.5rem',
              fontWeight: 'bold',
              cursor: 'pointer',
              transition: 'all 0.2s ease'
            }}
          >
            🎤 Start Listening
          </button>
        ) : isListening ? (
          <button
            className="btn"
            onClick={handleStopListening}
            style={{
              fontSize: '0.875rem',
              padding: '0.5rem 1rem',
              backgroundColor: '#ff4444',
              color: 'var(--paper-white)',
              border: '2px solid #ff6666',
              borderRadius: '0.5rem',
              fontWeight: 'bold',
              cursor: 'pointer',
              transition: 'all 0.2s ease'
            }}
          >
            🛑 Stop
          </button>
        ) : (
          <div style={{
            fontSize: '0.875rem',
            padding: '0.5rem 1rem',
            opacity: 0.7,
            fontStyle: 'italic',
            color: '#cccccc'
          }}>
            {isRestartingVoice ? '🔄 Restarting...' : 
             isSpeaking ? '🎙️ Generating...' :
             isPlaying ? '🔊 Speaking...' : 
             conversationActive ? '✅ Ready for next question' :
             '⏳ System Ready'}
          </div>
        )}
      </div>
      
      {/* Visual Listening Indicator */}
      {isListening && (
        <div style={{
          display: 'flex',
          justifyContent: 'center',
          gap: '4px',
          marginTop: '0.5rem'
        }}>
          {[1, 2, 3, 4, 5].map((i) => (
            <div
              key={i}
              style={{
                width: '3px',
                height: '16px',
                backgroundColor: 'var(--paper-white)',
                borderRadius: '2px',
                animation: `pulse 1.5s ease-in-out ${i * 0.2}s infinite`
              }}
            />
          ))}
        </div>
      )}

      {/* Processing Animation */}
      {(isSpeaking || isPlaying || isRestartingVoice) && (
        <div style={{
          display: 'flex',
          justifyContent: 'center',
          gap: '6px',
          marginTop: '0.5rem'
        }}>
          {[1, 2, 3].map((i) => (
            <div
              key={i}
              style={{
                width: '8px',
                height: '8px',
                backgroundColor: 'var(--paper-white)',
                borderRadius: '50%',
                animation: `bounce 1.4s ease-in-out ${i * 0.16}s infinite both`
              }}
            />
          ))}
        </div>
      )}
    </div>
  )
}

export default VoiceInterface 