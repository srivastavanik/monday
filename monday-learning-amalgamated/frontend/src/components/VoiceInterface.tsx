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

  // Smart button logic - prevent conflicts
  const canManuallyStart = !isListening && !isRestartingVoice && !isSpeaking && !isPlaying
  const shouldShowStartButton = canManuallyStart && conversationActive
  
  const handleStartListening = () => {
    // Double-check conditions before attempting start
    if (!canManuallyStart) {
      console.log('🚫 Manual start blocked - system busy')
      return
    }
    
    onUserInteraction?.()
    onStartListening()
  }

  const handleStopListening = () => {
    onUserInteraction?.()
    onStopListening()
  }

  // Status message logic
  const getStatusMessage = () => {
    if (isSpeaking || isPlaying) {
      return 'Monday is speaking...'
    }
    if (isRestartingVoice) {
      return 'Restarting voice recognition...'
    }
    if (isListening) {
      return conversationActive ? 'Listening...' : 'Listening for "Hey Monday..."'
    }
    if (conversationActive) {
      return 'In conversation - voice ready'
    }
    return 'Voice recognition ready'
  }

  return (
    <div style={{
      position: 'fixed',
      bottom: '2rem',
      left: '50%',
      transform: 'translateX(-50%)',
      padding: '1rem 2rem',
      backgroundColor: isListening ? 'var(--true-turquoise)' : 'rgba(9, 23, 23, 0.8)',
      color: 'var(--paper-white)',
      borderRadius: '2rem',
      border: '2px solid var(--true-turquoise)',
      minWidth: '300px',
      textAlign: 'center',
      zIndex: 1000,
      backdropFilter: 'blur(10px)'
    }}>
      <div style={{
        fontSize: '0.875rem',
        marginBottom: '0.5rem',
        opacity: 0.8
      }}>
        {getStatusMessage()}
      </div>
      
      {transcript && (
        <div style={{
          fontSize: '1rem',
          fontWeight: 500,
          padding: '0.5rem 0'
        }}>
          "{transcript}"
        </div>
      )}
      
      {conversationActive && (
        <div style={{
          fontSize: '0.75rem',
          color: 'var(--true-turquoise)',
          marginBottom: '0.5rem'
        }}>
          💬 You can speak without saying "Monday"
        </div>
      )}
      
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
              padding: '0.5rem 1rem'
            }}
          >
            Start Listening
          </button>
        ) : isListening ? (
          <button
            className="btn"
            onClick={handleStopListening}
            style={{
              fontSize: '0.875rem',
              padding: '0.5rem 1rem',
              backgroundColor: 'var(--paper-white)',
              color: 'var(--true-turquoise)'
            }}
          >
            Stop
          </button>
        ) : (
          <div style={{
            fontSize: '0.875rem',
            padding: '0.5rem 1rem',
            opacity: 0.7,
            fontStyle: 'italic'
          }}>
            {isRestartingVoice ? 'Restarting...' : 
             isSpeaking || isPlaying ? 'Speaking...' : 
             'Voice Ready'}
          </div>
        )}
      </div>
      
      {isListening && (
        <div style={{
          display: 'flex',
          justifyContent: 'center',
          gap: '4px',
          marginTop: '0.5rem'
        }}>
          {[1, 2, 3].map((i) => (
            <div
              key={i}
              style={{
                width: '4px',
                height: '20px',
                backgroundColor: 'var(--paper-white)',
                borderRadius: '2px',
                animation: `pulse 1.5s ease-in-out ${i * 0.3}s infinite`
              }}
            />
          ))}
        </div>
      )}
    </div>
  )
}

export default VoiceInterface 