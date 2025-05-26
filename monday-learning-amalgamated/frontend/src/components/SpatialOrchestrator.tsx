import { useRef } from 'react'
import { useMondayStore } from '../store/mondayStore'
import InformationPanel from './InformationPanel'

const SpatialOrchestrator = () => {
  const groupRef = useRef<any>()
  const {
    panels,
    activePanel,
    setActivePanel
  } = useMondayStore()

  // Handle panel selection
  const handlePanelSelect = (panelId: string) => {
    setActivePanel(panelId === activePanel ? null : panelId)
  }

  // Calculate static positions for panels
  const getPanelPosition = (index: number): [number, number, number] => {
    if (panels.length === 1) {
      // Single panel centered
      return [0, 1.6, -2]
    } else if (panels.length === 2) {
      // Two panels side by side
      return index === 0 ? [-1.5, 1.6, -2] : [1.5, 1.6, -2]
    } else {
      // Multiple panels in a horizontal line
      const spacing = 2.5
      const startX = -(panels.length - 1) * spacing / 2
      return [startX + index * spacing, 1.6, -2]
    }
  }

  return (
    <group ref={groupRef}>
      {/* Render information panels */}
      {panels.map((panel, index) => {
        const position = getPanelPosition(index)
        
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
    </group>
  )
}

export default SpatialOrchestrator 