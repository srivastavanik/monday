import { useRef, useEffect } from 'react'
import { useFrame } from '@react-three/fiber'
import { useMondayStore } from '../store/mondayStore'
import InformationPanel from './InformationPanel'

const SpatialOrchestrator = () => {
  const groupRef = useRef<any>()
  const {
    panels,
    activePanel,
    reasoningChain,
    researchWeb,
    focusMode,
    spatialLayout,
    setActivePanel
  } = useMondayStore()

  // Debug logging
  useEffect(() => {
    console.log('SpatialOrchestrator: Panels updated:', panels.length, panels)
  }, [panels])

  // Automatic spatial layout management
  useFrame(() => {
    if (!groupRef.current) return

    // Rotate the entire constellation slowly
    groupRef.current.rotation.y += 0.001

    // Don't mutate the original panel positions - create local positions instead
    // The position mutation was causing issues with React's rendering
  })

  // Handle panel selection
  const handlePanelSelect = (panelId: string) => {
    setActivePanel(panelId === activePanel ? null : panelId)
  }

  // Calculate positions for panels without mutating the original data
  const getPanelPosition = (panel: any, index: number): [number, number, number] => {
    if (spatialLayout === 'focus' && panel.id !== activePanel) {
      // Move non-active panels away in focus mode
      const distance = 5
      const angle = (index / panels.length) * Math.PI * 2
      return [Math.cos(angle) * distance, panel.position[1], Math.sin(angle) * distance]
    } else if (spatialLayout === 'research') {
      // Spread panels in research constellation
      const radius = 3
      const layers = Math.ceil(panels.length / 6)
      const layer = Math.floor(index / 6)
      const angleStep = (Math.PI * 2) / Math.min(6, panels.length - layer * 6)
      const angle = (index % 6) * angleStep
      
      return [
        Math.cos(angle) * (radius + layer * 1.5),
        1.6 + layer * 0.5,
        Math.sin(angle) * (radius + layer * 1.5)
      ]
    } else {
      // Default semicircle layout
      const angle = (index / Math.max(panels.length - 1, 1)) * Math.PI - Math.PI / 2
      const radius = 2.5
      
      return [
        Math.cos(angle) * radius,
        1.6 + Math.sin(index * 0.5) * 0.3,
        Math.sin(angle) * radius * 0.5 - 1
      ]
    }
  }

  console.log('SpatialOrchestrator: Rendering with', panels.length, 'panels')

  return (
    <group ref={groupRef}>
      {/* Render information panels */}
      {panels.map((panel, index) => {
        const position = getPanelPosition(panel, index)
        console.log(`SpatialOrchestrator: Rendering panel ${index + 1}:`, panel.title, 'at position:', position)
        
        return (
          <InformationPanel
            key={panel.id}
            panel={{
              ...panel,
              position: position
            }}
            isActive={panel.id === activePanel}
            onSelect={() => handlePanelSelect(panel.id)}
          />
        )
      })}

      {/* Render reasoning chain visualization */}
      {reasoningChain.map((step, index) => (
        <group key={step.id} position={step.position}>
          {/* Reasoning step sphere */}
          <mesh position={[0, 0, 0]}>
            <sphereGeometry args={[0.1]} />
            <meshStandardMaterial 
              color="#20808D" 
              transparent 
              opacity={step.confidence}
            />
          </mesh>
          
          {/* Connection lines to next step */}
          {index < reasoningChain.length - 1 && (
            <line>
              <bufferGeometry>
                <bufferAttribute
                  attach="attributes-position"
                  count={2}
                  array={new Float32Array([
                    ...step.position,
                    ...reasoningChain[index + 1].position
                  ])}
                  itemSize={3}
                />
              </bufferGeometry>
              <lineBasicMaterial color="#20808D" opacity={0.5} transparent />
            </line>
          )}
        </group>
      ))}

      {/* Render research web nodes */}
      {researchWeb.map((node) => (
        <group key={node.id} position={node.position}>
          <mesh>
            <boxGeometry args={[0.15, 0.15, 0.15]} />
            <meshStandardMaterial 
              color={node.type === 'source' ? '#FBFAF4' : '#20808D'}
              transparent
              opacity={0.8}
            />
          </mesh>
          
          {/* Render connections between nodes */}
          {node.connections.map((connectionId) => {
            const targetNode = researchWeb.find(n => n.id === connectionId)
            if (!targetNode) return null
            
            return (
              <line key={connectionId}>
                <bufferGeometry>
                  <bufferAttribute
                    attach="attributes-position"
                    count={2}
                    array={new Float32Array([
                      ...node.position,
                      ...targetNode.position
                    ])}
                    itemSize={3}
                  />
                </bufferGeometry>
                <lineBasicMaterial color="#20808D" opacity={0.3} transparent />
              </line>
            )
          })}
        </group>
      ))}

      {/* Ambient particles for atmosphere */}
      <group>
        {Array.from({ length: 20 }, (_, i) => (
          <mesh
            key={i}
            position={[
              (Math.random() - 0.5) * 20,
              Math.random() * 5,
              (Math.random() - 0.5) * 20
            ]}
          >
            <sphereGeometry args={[0.01]} />
            <meshBasicMaterial 
              color="#20808D" 
              transparent 
              opacity={0.2}
            />
          </mesh>
        ))}
      </group>
    </group>
  )
}

export default SpatialOrchestrator 