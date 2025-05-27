import React, { useState, useEffect, useRef } from 'react'

interface ProgressivePanelProps {
  id: string
  title: string
  isVisible: boolean
  onClose: () => void
}

interface ProgressUpdate {
  type: string
  message: string
  progress: number
  reasoning?: string
  sources?: string[]
}

const ProgressivePanel: React.FC<ProgressivePanelProps> = ({
  id,
  title,
  isVisible,
  onClose
}) => {
  const [progressData, setProgressData] = useState<ProgressUpdate | null>(null)
  const [displayContent, setDisplayContent] = useState<string>('')
  const [isComplete, setIsComplete] = useState(false)
  const socketRef = useRef<any>(null)

  useEffect(() => {
    // Access the global socket
    const socket = (window as any).socket
    if (!socket) {
      console.warn('🔄 ProgressivePanel: No global socket available')
      return
    }

    socketRef.current = socket
    console.log('🔄 ProgressivePanel: Connected to socket for real-time updates')

    const handleReasoningProgress = (data: any) => {
      console.log('🔄 ProgressivePanel: Received reasoning progress:', data)
      
      if (data.update) {
        setProgressData(data.update)
        
        // Update display content with paragraph-level rotation
        if (data.update.reasoning && data.update.reasoning.length > 50) {
          let cleanContent = data.update.reasoning
          
          // Clean up thinking tags
          if (cleanContent.startsWith('<think>')) {
            const thinkMatch = cleanContent.match(/<think>([\s\S]*?)<\/think>/);
            if (thinkMatch) {
              cleanContent = thinkMatch[1].trim()
            } else {
              cleanContent = cleanContent.replace(/^<think>\s*/g, '').trim()
            }
          }
          
          // Show latest paragraphs for progressive display
          const paragraphs = cleanContent.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 30)
          
          if (paragraphs.length > 0) {
            let displayText = ''
            if (paragraphs.length === 1) {
              displayText = paragraphs[0].trim()
            } else if (paragraphs.length >= 2) {
              const lastTwo = paragraphs.slice(-2).join('\n\n').trim()
              displayText = lastTwo.length > 1000 ? paragraphs[paragraphs.length - 1].trim() : lastTwo
            }
            
            // Add continuation indicator if there are more paragraphs
            if (paragraphs.length > 2) {
              displayText = `[Analysis in progress...]\n\n${displayText}`
            }
            
            setDisplayContent(displayText)
          }
        }
        
        // Mark as complete when progress reaches 100%
        if (data.update.progress === 100 || data.update.type === 'complete') {
          setIsComplete(true)
          console.log('🔄 ProgressivePanel: Research completed')
        }
      }
    }

    const handleResearchProgress = (data: any) => {
      console.log('🔄 ProgressivePanel: Received research progress:', data)
      handleReasoningProgress(data) // Use same handler for both
    }

    // Listen for progress updates
    socket.on('reasoning_progress', handleReasoningProgress)
    socket.on('research_progress', handleResearchProgress)

    return () => {
      if (socketRef.current) {
        socketRef.current.off('reasoning_progress', handleReasoningProgress)
        socketRef.current.off('research_progress', handleResearchProgress)
      }
    }
  }, [])

  if (!isVisible) return null

  return (
    <div style={{
      position: 'fixed',
      top: '50%',
      right: '20px',
      transform: 'translateY(-50%)',
      width: '400px',
      maxHeight: '70vh',
      backgroundColor: 'var(--paper-white)',
      border: '2px solid var(--true-turquoise)',
      borderRadius: '8px',
      boxShadow: '0 8px 32px rgba(32, 128, 141, 0.3)',
      zIndex: 1000,
      overflow: 'hidden',
      display: 'flex',
      flexDirection: 'column'
    }}>
      {/* Header */}
      <div style={{
        padding: '1rem',
        backgroundColor: 'var(--true-turquoise)',
        color: 'var(--paper-white)',
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center'
      }}>
        <h3 style={{ margin: 0, fontSize: '1rem' }}>{title}</h3>
        <button
          onClick={onClose}
          style={{
            background: 'none',
            border: 'none',
            color: 'var(--paper-white)',
            fontSize: '1.2rem',
            cursor: 'pointer',
            padding: '0',
            width: '24px',
            height: '24px',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center'
          }}
        >
          ×
        </button>
      </div>

      {/* Progress Bar */}
      {progressData && !isComplete && (
        <div style={{
          padding: '0.5rem 1rem',
          backgroundColor: 'rgba(32, 128, 141, 0.1)',
          borderBottom: '1px solid rgba(32, 128, 141, 0.2)'
        }}>
          <div style={{
            display: 'flex',
            justifyContent: 'space-between',
            alignItems: 'center',
            marginBottom: '0.5rem'
          }}>
            <span style={{ fontSize: '0.875rem', color: 'var(--true-turquoise)' }}>
              {progressData.message}
            </span>
            <span style={{ fontSize: '0.875rem', fontWeight: 'bold', color: 'var(--true-turquoise)' }}>
              {progressData.progress}%
            </span>
          </div>
          <div style={{
            width: '100%',
            height: '4px',
            backgroundColor: 'rgba(32, 128, 141, 0.2)',
            borderRadius: '2px',
            overflow: 'hidden'
          }}>
            <div style={{
              width: `${progressData.progress}%`,
              height: '100%',
              backgroundColor: 'var(--true-turquoise)',
              transition: 'width 0.3s ease'
            }} />
          </div>
        </div>
      )}

      {/* Content */}
      <div style={{
        flex: 1,
        padding: '1rem',
        overflow: 'auto',
        fontSize: '0.875rem',
        lineHeight: '1.5',
        color: 'var(--offblack)'
      }}>
        {isComplete && (
          <div style={{
            padding: '0.5rem',
            backgroundColor: 'rgba(34, 197, 94, 0.1)',
            border: '1px solid rgba(34, 197, 94, 0.3)',
            borderRadius: '4px',
            marginBottom: '1rem',
            color: 'rgb(34, 197, 94)',
            fontSize: '0.875rem',
            fontWeight: 'bold'
          }}>
            ✅ Analysis Complete
          </div>
        )}
        
        {displayContent ? (
          <div style={{ whiteSpace: 'pre-wrap' }}>
            {displayContent}
          </div>
        ) : (
          <div style={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            height: '100px',
            color: 'var(--true-turquoise)'
          }}>
            <div style={{
              width: '20px',
              height: '20px',
              border: '2px solid var(--true-turquoise)',
              borderTop: '2px solid transparent',
              borderRadius: '50%',
              animation: 'spin 1s linear infinite',
              marginRight: '0.5rem'
            }} />
            Initializing research...
          </div>
        )}
        
        {progressData?.sources && progressData.sources.length > 0 && (
          <div style={{
            marginTop: '1rem',
            padding: '0.5rem',
            backgroundColor: 'rgba(32, 128, 141, 0.1)',
            borderRadius: '4px'
          }}>
            <div style={{ fontWeight: 'bold', marginBottom: '0.5rem', color: 'var(--true-turquoise)' }}>
              Sources:
            </div>
            <ul style={{ margin: 0, paddingLeft: '1rem' }}>
              {progressData.sources.map((source, index) => (
                <li key={index} style={{ marginBottom: '0.25rem' }}>{source}</li>
              ))}
            </ul>
          </div>
        )}
      </div>
    </div>
  )
}

export default ProgressivePanel 