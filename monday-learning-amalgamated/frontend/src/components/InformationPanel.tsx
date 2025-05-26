import React, { useRef, useMemo } from 'react'
import { Text, Plane } from '@react-three/drei'
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

  // Panel dimensions based on content
  const dimensions = useMemo(() => {
    const baseWidth = 2.5
    const baseHeight = 1.8
    const contentLength = panel.content.length
    
    // Adjust size based on content
    const width = Math.min(baseWidth + (contentLength / 1000) * 0.5, 4.0)
    const height = Math.min(baseHeight + (contentLength / 500) * 0.3, 3.5)
    
    return { width, height }
  }, [panel.content])

  // Color scheme based on panel type and app theming
  const colors = useMemo(() => {
    const baseColors = {
      content: '#FBFAF4',      // Paper White
      reasoning: '#20808D',     // True Turquoise
      research: '#091717',      // Offblack
      visualization: '#20808D', // True Turquoise
      citations: '#FBFAF4'      // Paper White
    }
    
    return {
      background: baseColors[panel.type] || baseColors.content,
      text: panel.type === 'research' ? '#FBFAF4' : '#091717',
      border: '#20808D',
      glow: isActive ? '#20808D' : 'transparent'
    }
  }, [panel.type, isActive])

  // Format content for display - show "Monday:" prefix
  const formattedContent = useMemo(() => {
    let content = panel.content
    if (content.length > 500) {
      content = content.substring(0, 497) + '...'
    }
    
    // Add "Monday:" prefix if not already present
    if (!content.toLowerCase().startsWith('monday:')) {
      content = `Monday: ${content}`
    }
    
    return content
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
        position={[0, dimensions.height/2 - 0.15, 0]}
        fontSize={0.1}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        maxWidth={dimensions.width - 0.2}
      >
        {panel.title}
      </Text>
      
      {/* Panel content */}
      <Text
        position={[0, -0.1, 0]}
        fontSize={0.06}
        color={colors.text}
        anchorX="center"
        anchorY="middle"
        maxWidth={dimensions.width - 0.3}
        lineHeight={1.3}
      >
        {formattedContent}
      </Text>
      
      {/* Citations indicator */}
      {panel.citations && panel.citations.length > 0 && (
        <group position={[dimensions.width/2 - 0.15, -dimensions.height/2 + 0.1, 0.01]}>
          <Plane args={[0.25, 0.08]}>
            <meshStandardMaterial 
              color={colors.border}
              transparent
              opacity={0.8}
            />
          </Plane>
          <Text
            position={[0, 0, 0.01]}
            fontSize={0.04}
            color="#FBFAF4"
            anchorX="center"
            anchorY="middle"
          >
            {panel.citations.length} sources
          </Text>
        </group>
      )}
      
      {/* Type indicator */}
      <group position={[-dimensions.width/2 + 0.15, dimensions.height/2 - 0.08, 0.01]}>
        <Plane args={[0.25, 0.08]}>
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
          {panel.type.toUpperCase()}
        </Text>
      </group>
    </group>
  )
}

export default InformationPanel 