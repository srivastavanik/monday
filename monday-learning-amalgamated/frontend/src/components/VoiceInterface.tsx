interface VoiceInterfaceProps {
  isListening: boolean
  transcript: string
  onStartListening: () => void
  onStopListening: () => void
}

const VoiceInterface = ({ 
  isListening, 
  transcript, 
  onStartListening, 
  onStopListening 
}: VoiceInterfaceProps) => {
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
        {isListening ? 'Listening for "Monday..."' : 'Voice recognition ready'}
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
      
      <div style={{
        display: 'flex',
        justifyContent: 'center',
        gap: '1rem',
        marginTop: '0.5rem'
      }}>
        {!isListening ? (
          <button
            className="btn"
            onClick={onStartListening}
            style={{
              fontSize: '0.875rem',
              padding: '0.5rem 1rem'
            }}
          >
            Start Listening
          </button>
        ) : (
          <button
            className="btn"
            onClick={onStopListening}
            style={{
              fontSize: '0.875rem',
              padding: '0.5rem 1rem',
              backgroundColor: 'var(--paper-white)',
              color: 'var(--true-turquoise)'
            }}
          >
            Stop
          </button>
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