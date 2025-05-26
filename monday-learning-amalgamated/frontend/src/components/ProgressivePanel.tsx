import React, { useState, useEffect, useRef } from 'react'
import { Text, Plane } from '@react-three/drei'

interface ProgressUpdate {
  type: 'thinking' | 'researching' | 'searching' | 'analyzing' | 'synthesizing' | 'complete'
  message: string
  progress: number
  sources?: string[]
  reasoning?: string
  metadata?: any
}

interface ProgressivePanelProps {
  id: string
  title: string
  position: [number, number, number]
  rotation: [number, number, number]
  isActive: boolean
  onSelect?: () => void
  type: 'progressive_reasoning' | 'progressive_research'
  initialContent?: string
}

const ProgressivePanel: React.FC<ProgressivePanelProps> = ({
  id,
  title,
  position,
  rotation,
  isActive,
  onSelect,
  type,
  initialContent = ''
}) => {
  const [currentContent, setCurrentContent] = useState(initialContent)
  const [progress, setProgress] = useState(0)
  const [currentMessage, setCurrentMessage] = useState('Starting analysis...')
  const [sources, setSources] = useState<string[]>([])
  const [isComplete, setIsComplete] = useState(false)
  const meshRef = useRef<any>()

  // Update content when initialContent changes
  useEffect(() => {
    if (initialContent && initialContent.length > 0) {
      console.log('📊 ProgressivePanel: Setting initial content:', {
        contentLength: initialContent.length,
        contentPreview: initialContent.substring(0, 100)
      })
      setCurrentContent(initialContent)
      setProgress(10) // Show some initial progress
      setCurrentMessage('Processing your request...')
    }
  }, [initialContent])

  // Listen for progress updates from the backend
  useEffect(() => {
    // Get socket from the global window object (set by App.tsx)
    const socket = (window as any).socket
    if (!socket) {
      console.warn('📊 ProgressivePanel: No socket available for progress updates')
      return
    }

    console.log('📊 ProgressivePanel: Setting up progress listeners for', type, 'with socket ID:', socket.id)

    const handleProgressUpdate = (data: any) => {
      console.log('📊 ProgressivePanel: Received progress update:', {
        type: data.type,
        updateType: data.update?.type,
        progress: data.update?.progress,
        message: data.update?.message,
        panelType: type,
        hasReasoning: !!data.update?.reasoning,
        socketId: socket.id
      })
      
      if (data.type === 'progress_update') {
        const update: ProgressUpdate = data.update
        
        console.log('📊 ProgressivePanel: Processing update:', {
          type: update.type,
          progress: update.progress,
          messageLength: update.message?.length || 0,
          reasoningLength: update.reasoning?.length || 0
        })
        
        setProgress(update.progress || 0)
        setCurrentMessage(update.message || 'Processing...')
        
        if (update.sources) {
          setSources(update.sources)
        }
        
        // Update content with reasoning if available
        if (update.reasoning && update.reasoning.length > 0) {
          console.log('📊 ProgressivePanel: Updating content with reasoning:', update.reasoning.substring(0, 100))
          setCurrentContent(update.reasoning)
        }
        
        if (update.type === 'complete') {
          setIsComplete(true)
          if (update.reasoning) {
            setCurrentContent(update.reasoning)
          }
          console.log('📊 ProgressivePanel: Process completed, final content length:', update.reasoning?.length || 0)
        }
      }
    }

    // Listen for both reasoning and research progress
    socket.on('reasoning_progress', handleProgressUpdate)
    socket.on('research_progress', handleProgressUpdate)

    console.log('📊 ProgressivePanel: Progress listeners set up successfully for socket:', socket.id)

    return () => {
      console.log('📊 ProgressivePanel: Cleaning up progress listeners for socket:', socket.id)
      socket.off('reasoning_progress', handleProgressUpdate)
      socket.off('research_progress', handleProgressUpdate)
    }
  }, [type])

  // Debug: Log when content changes
  useEffect(() => {
    console.log('📊 ProgressivePanel: Content updated:', {
      id,
      type,
      contentLength: currentContent.length,
      progress,
      isComplete,
      currentMessage
    })
  }, [currentContent, progress, isComplete, currentMessage, id, type])

  // Panel dimensions
  const width = 3.0
  const height = 2.5

  // Color scheme based on type
  const colors = {
    background: type === 'progressive_reasoning' ? '#20808D' : '#091717',
    text: type === 'progressive_reasoning' ? '#FBFAF4' : '#FBFAF4',
    border: '#20808D',
    progress: '#FBFAF4',
    glow: isActive ? '#20808D' : 'transparent'
  }

  // Format content for display
  const displayContent = currentContent.length > 800 
    ? currentContent.substring(0, 797) + '...'
    : currentContent

  // Progress bar width
  const progressBarWidth = (width - 0.4) * (progress / 100)

  return (
    <group 
      ref={meshRef}
      position={position}
      rotation={rotation}
      onClick={onSelect}
    >
      {/* Panel background */}
      <Plane args={[width, height]} position={[0, 0, -0.01]}>
        <meshStandardMaterial 
          color={colors.background}
          transparent
          opacity={0.9}
        />
      </Plane>
      
      {/* Border */}
      <Plane args={[width + 0.05, height + 0.05]} position={[0, 0, -0.02]}>
        <meshStandardMaterial 
          color={colors.border}
          transparent
          opacity={0.8}
        />
      </Plane>
      
      {/* Glow effect for active panels */}
      {isActive && (
        <Plane args={[width + 0.2, height + 0.2]} position={[0, 0, -0.03]}>
          <meshStandardMaterial 
            color={colors.glow}
            transparent
            opacity={0.3}
          />
        </Plane>
      )}
      
      {/* Panel title */}
      <Text
        position={[0, height/2 - 0.15, 0]}
        fontSize={0.1}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        maxWidth={width - 0.2}
      >
        {title}
      </Text>
      
      {/* Progress bar background */}
      <Plane 
        args={[width - 0.4, 0.08]} 
        position={[0, height/2 - 0.35, 0.001]}
      >
        <meshStandardMaterial 
          color="#333333"
          transparent
          opacity={0.6}
        />
      </Plane>
      
      {/* Progress bar fill */}
      <Plane 
        args={[progressBarWidth, 0.06]} 
        position={[-(width - 0.4)/2 + progressBarWidth/2, height/2 - 0.35, 0.002]}
      >
        <meshStandardMaterial 
          color={colors.progress}
          transparent
          opacity={0.9}
        />
      </Plane>
      
      {/* Progress percentage */}
      <Text
        position={[width/2 - 0.2, height/2 - 0.35, 0.003]}
        fontSize={0.05}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
      >
        {Math.round(progress)}%
      </Text>
      
      {/* Current status message */}
      <Text
        position={[0, height/2 - 0.55, 0]}
        fontSize={0.06}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        maxWidth={width - 0.3}
      >
        {currentMessage}
      </Text>
      
      {/* Main content area */}
      <Text
        position={[0, -0.1, 0]}
        fontSize={0.05}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        maxWidth={width - 0.3}
        lineHeight={1.2}
      >
        {displayContent}
      </Text>
      
      {/* Sources indicator */}
      {sources.length > 0 && (
        <group position={[0, -height/2 + 0.2, 0.01]}>
          <Text
            position={[0, 0, 0]}
            fontSize={0.04}
            color={colors.text}
            anchorX="center"
            anchorY="middle"
            maxWidth={width - 0.2}
          >
            Sources: {sources.join(', ')}
          </Text>
        </group>
      )}
      
      {/* Completion indicator */}
      {isComplete && (
        <group position={[width/2 - 0.15, height/2 - 0.08, 0.01]}>
          <Plane args={[0.25, 0.08]}>
            <meshStandardMaterial 
              color="#00ff00"
              transparent
              opacity={0.8}
            />
          </Plane>
          <Text
            position={[0, 0, 0.01]}
            fontSize={0.04}
            color="#000000"
            anchorX="center"
            anchorY="middle"
          >
            COMPLETE
          </Text>
        </group>
      )}
      
      {/* Type indicator */}
      <group position={[-width/2 + 0.2, height/2 - 0.08, 0.01]}>
        <Plane args={[0.35, 0.08]}>
          <meshStandardMaterial 
            color={colors.border}
            transparent
            opacity={0.6}
          />
        </Plane>
        <Text
          position={[0, 0, 0.01]}
          fontSize={0.035}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {type === 'progressive_reasoning' ? 'REASONING' : 'RESEARCH'}
        </Text>
      </group>
    </group>
  )
}

export default ProgressivePanel 