import React, { memo } from 'react'
import { useMondayStore } from '../store/mondayStore'
import InformationPanel from './InformationPanel'

const SpatialOrchestrator = memo(() => {
  const { panels } = useMondayStore()

  // Only render the first two panels in fixed positions
  const mainPanel = panels[0]
  const reasoningPanel = panels[1]

  return (
    <group>
      {/* Main panel on the left */}
      {mainPanel && (
        <InformationPanel
          key={mainPanel.id}
          panel={{
            ...mainPanel,
            position: [-1.5, 1.6, -2] as [number, number, number]
          }}
          isActive={true}
          onSelect={() => {}}
        />
      )}
      
      {/* Reasoning panel on the right */}
      {reasoningPanel && (
        <InformationPanel
          key={reasoningPanel.id}
          panel={{
            ...reasoningPanel,
            position: [1.5, 1.6, -2] as [number, number, number]
          }}
          isActive={false}
          onSelect={() => {}}
        />
      )}
    </group>
  )
})

SpatialOrchestrator.displayName = 'SpatialOrchestrator'

export default SpatialOrchestrator 