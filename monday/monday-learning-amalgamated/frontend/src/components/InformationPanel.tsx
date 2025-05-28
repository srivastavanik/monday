import React, { useRef, useMemo } from 'react'
import { useFrame } from '@react-three/fiber'
import { Text, Plane, Html } from '@react-three/drei'
import { MondayPanel } from '../store/mondayStore'

interface InformationPanelProps {
  panel: MondayPanel
  isActive: boolean
  onSelect?: () => void
}

const InformationPanel: React.FC<InformationPanelProps> = ({ 
  panel, 
  isActive, 
  onSelect 
}) => {
  const meshRef = useRef<any>()
  const textRef = useRef<any>()
  
  // Animate panel based on state
  useFrame((state) => {
    if (meshRef.current) {
      // Gentle floating animation
      meshRef.current.position.y = panel.position[1] + Math.sin(state.clock.elapsedTime * 0.5) * 0.05
      
      // Glow effect for active panels
      if (isActive) {
        meshRef.current.scale.setScalar(1 + Math.sin(state.clock.elapsedTime * 2) * 0.02)
      }
    }
  })

  // Panel dimensions based on content
  const dimensions = useMemo(() => {
    const baseWidth = 1.5
    const baseHeight = 1.0
    const contentLength = panel.content.length
    
    // Adjust size based on content
    const width = Math.min(baseWidth + (contentLength / 1000) * 0.5, 3.0)
    const height = Math.min(baseHeight + (contentLength / 500) * 0.3, 2.0)
    
    return { width, height }
  }, [panel.content])

  // Color scheme based on panel type and Perplexity branding
  const colors = useMemo(() => {
    const baseColors = {
      content: '#FBFAF4',      // Paper White
      reasoning: '#20808D',     // True Turquoise
      research: '#091717',      // Offblack
      visualization: '#20808D'  // True Turquoise
    }
    
    return {
      background: baseColors[panel.type] || baseColors.content,
      text: panel.type === 'research' ? '#FBFAF4' : '#091717',
      border: '#20808D',
      glow: isActive ? '#20808D' : 'transparent'
    }
  }, [panel.type, isActive])

  // Format content for display
  const formattedContent = useMemo(() => {
    if (panel.content.length > 300) {
      return panel.content.substring(0, 297) + '...'
    }
    return panel.content
  }, [panel.content])

  return (
    <group 
      ref={meshRef}
      position={panel.position}
      rotation={panel.rotation}
      onClick={onSelect}
    >
      {/* Panel background */}
      <Plane args={[dimensions.width, dimensions.height]} position={[0, 0, -0.01]}>
        <meshStandardMaterial 
          color={colors.background}
          transparent
          opacity={panel.opacity}
        />
      </Plane>
      
      {/* Border */}
      <Plane args={[dimensions.width + 0.05, dimensions.height + 0.05]} position={[0, 0, -0.02]}>
        <meshStandardMaterial 
          color={colors.border}
          transparent
          opacity={panel.opacity * 0.8}
        />
      </Plane>
      
      {/* Glow effect for active panels */}
      {isActive && (
        <Plane args={[dimensions.width + 0.2, dimensions.height + 0.2]} position={[0, 0, -0.03]}>
          <meshStandardMaterial 
            color={colors.glow}
            transparent
            opacity={0.3}
          />
        </Plane>
      )}
      
      {/* Panel title */}
      <Text
        position={[0, dimensions.height/2 - 0.1, 0]}
        fontSize={0.08}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        font="/fonts/FK-Display-Medium.woff"
        maxWidth={dimensions.width - 0.2}
      >
        {panel.title}
      </Text>
      
      {/* Panel content */}
      <Text
        ref={textRef}
        position={[0, 0, 0]}
        fontSize={0.05}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        font="/fonts/FK-Grotesk-Neue-Regular.woff"
        maxWidth={dimensions.width - 0.2}
        lineHeight={1.2}
      >
        {formattedContent}
      </Text>
      
      {/* Citations indicator */}
      {panel.citations && panel.citations.length > 0 && (
        <group position={[dimensions.width/2 - 0.1, -dimensions.height/2 + 0.1, 0.01]}>
          <Plane args={[0.15, 0.05]}>
            <meshStandardMaterial 
              color={colors.border}
              transparent
              opacity={0.8}
            />
          </Plane>
          <Text
            position={[0, 0, 0.01]}
            fontSize={0.03}
            color="#FBFAF4"
            anchorX="center"
            anchorY="middle"
          >
            {panel.citations.length} refs
          </Text>
        </group>
      )}
      
      {/* Type indicator */}
      <group position={[-dimensions.width/2 + 0.1, dimensions.height/2 - 0.05, 0.01]}>
        <Plane args={[0.2, 0.05]}>
          <meshStandardMaterial 
            color={colors.border}
            transparent
            opacity={0.6}
          />
        </Plane>
        <Text
          position={[0, 0, 0.01]}
          fontSize={0.025}
          color="#FBFAF4"
          anchorX="center"
          anchorY="middle"
        >
          {panel.type.toUpperCase()}
        </Text>
      </group>
      
      {/* Interactive hover effect */}
      <Plane 
        args={[dimensions.width, dimensions.height]} 
        position={[0, 0, 0.01]}
        visible={false}
      >
        <meshStandardMaterial transparent opacity={0} />
      </Plane>
    </group>
  )
}

export default InformationPanel 