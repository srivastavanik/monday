import React, { useRef, useMemo, memo, useCallback } from 'react'
import { useMondayStore } from '../store/mondayStore'
import InformationPanel from './InformationPanel'
import ProgressivePanel from './ProgressivePanel'

const SpatialOrchestrator = memo(() => {
  const groupRef = useRef<any>()
  const {
    panels,
    activePanel,
    setActivePanel
  } = useMondayStore()

  // Memoize panel selection handler to prevent re-renders
  const handlePanelSelect = useCallback((panelId: string) => {
    setActivePanel(panelId === activePanel ? null : panelId)
  }, [activePanel, setActivePanel])

  // Memoize panel positions to prevent infinite re-renders
  const panelPositions = useMemo(() => {
    return panels.map((panel, index) => {
      if (panels.length === 1) {
        // Single panel centered
        return [0, 1.6, -2] as [number, number, number]
      } else if (panels.length === 2) {
        // Two panels side by side - special handling for progressive panels
        if (panel.type === 'progressive_reasoning' || panel.type === 'progressive_research') {
          // Progressive panel on the right, closer and more visible
          return [2.2, 1.6, -1.5] as [number, number, number]
        } else {
          // Main panel on the left
          return [-1.2, 1.6, -2] as [number, number, number]
        }
      } else {
        // Multiple panels in a horizontal line
        const spacing = 2.5
        const startX = -(panels.length - 1) * spacing / 2
        return [startX + index * spacing, 1.6, -2] as [number, number, number]
      }
    })
  }, [panels.length, panels])

  // Memoize rendered panels to prevent unnecessary re-renders
  const renderedPanels = useMemo(() => {
    console.log('🎯 SpatialOrchestrator: Rendering panels, count:', panels.length)
    
    return panels.map((panel, index) => {
      const position = panelPositions[index]
      
      // Check if this is a progressive panel
      if (panel.type === 'progressive_reasoning' || panel.type === 'progressive_research') {
        return (
          <ProgressivePanel
            key={panel.id}
            id={panel.id}
            title={panel.title}
            position={position}
            rotation={panel.rotation}
            isActive={panel.id === activePanel}
            onSelect={() => handlePanelSelect(panel.id)}
            type={panel.type as 'progressive_reasoning' | 'progressive_research'}
            initialContent={panel.content || panel.fullContent || ''}
          />
        )
      }
      
      // Regular information panel
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
    })
  }, [panels, panelPositions, activePanel, handlePanelSelect])

  return (
    <group ref={groupRef}>
      {renderedPanels}
    </group>
  )
})

SpatialOrchestrator.displayName = 'SpatialOrchestrator'

export default SpatialOrchestrator 