import React from 'react'
import { SystemState, SystemStatus } from '../controllers/VoiceSystemController'

interface DiagnosticOverlayProps {
  systemState: SystemState
  systemStatus: SystemStatus
  transcript: string
  conversationActive: boolean
  error: string | null
  onEmergencyReset: () => void
  onForceStart: () => void
  onForceStop: () => void
}

const DiagnosticOverlay: React.FC<DiagnosticOverlayProps> = ({
  systemState,
  systemStatus,
  transcript,
  conversationActive,
  error,
  onEmergencyReset,
  onForceStart,
  onForceStop
}) => {
  const timeSinceLastTransition = Date.now() - systemStatus.lastTransition

  const getStateColor = (state: SystemState): string => {
    switch (state) {
      case SystemState.IDLE: return '#666666'
      case SystemState.WAITING_FOR_ACTIVATION: return '#00ff88'
      case SystemState.PROCESSING_COMMAND: return '#ffaa00'
      case SystemState.PLAYING_TTS: return '#ff6600'
      case SystemState.ACTIVE_LISTENING: return '#00ff00'
      case SystemState.ERROR: return '#ff0000'
      case SystemState.RESETTING: return '#ffff00'
      default: return '#ffffff'
    }
  }

  const getStatusIcon = (active: boolean): string => {
    return active ? '✅' : '❌'
  }

  // Detect potential issues
  const isStuck = timeSinceLastTransition > 10000
  const hasErrors = systemStatus.errorCount > 0
  const isInErrorState = systemState === SystemState.ERROR

  return (
    <div style={{
      position: 'fixed',
      top: '1rem',
      left: '1rem',
      padding: '1rem',
      backgroundColor: 'rgba(9, 23, 23, 0.95)',
      color: 'var(--paper-white)',
      fontSize: '0.75rem',
      borderRadius: '0.5rem',
      fontFamily: 'monospace',
      zIndex: 2000,
      border: `2px solid ${getStateColor(systemState)}`,
      maxWidth: '400px',
      minWidth: '350px'
    }}>
      {/* Header */}
      <div style={{
        fontWeight: 'bold',
        color: '#00ff88',
        marginBottom: '0.5rem',
        fontSize: '0.8rem',
        display: 'flex',
        justifyContent: 'space-between',
        alignItems: 'center'
      }}>
        <span>🔬 VOICE SYSTEM DIAGNOSTICS</span>
        <span style={{ 
          color: getStateColor(systemState),
          fontSize: '0.7rem',
          fontWeight: 'normal'
        }}>
          {systemState}
        </span>
      </div>

      {/* System State */}
      <div style={{ marginBottom: '0.5rem' }}>
        <div style={{ 
          color: getStateColor(systemState),
          fontWeight: 'bold',
          padding: '0.2rem 0.5rem',
          backgroundColor: 'rgba(255,255,255,0.1)',
          borderRadius: '0.25rem',
          marginBottom: '0.3rem'
        }}>
          State: {systemState}
        </div>
        
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0.2rem' }}>
          <div>Recognition: {getStatusIcon(systemStatus.recognitionActive)}</div>
          <div>TTS Gen: {getStatusIcon(systemStatus.ttsGenerating)}</div>
          <div>TTS Play: {getStatusIcon(systemStatus.ttsPlaying)}</div>
          <div>Conversation: {getStatusIcon(conversationActive)}</div>
        </div>
      </div>

      {/* Timing Information */}
      <div style={{ marginBottom: '0.5rem' }}>
        <div style={{ color: isStuck ? '#ff0000' : '#66ff66' }}>
          Last Transition: {Math.round(timeSinceLastTransition / 1000)}s ago
          {isStuck && ' ⚠️ STUCK!'}
        </div>
        <div style={{ color: hasErrors ? '#ff6666' : '#66ff66' }}>
          Error Count: {systemStatus.errorCount}
        </div>
      </div>

      {/* Current Activity */}
      {transcript && (
        <div style={{
          marginBottom: '0.5rem',
          padding: '0.3rem',
          backgroundColor: 'rgba(102,255,255,0.1)',
          borderRadius: '0.25rem',
          border: '1px solid #66ffff'
        }}>
          <div style={{ color: '#66ffff', fontSize: '0.7rem', marginBottom: '0.2rem' }}>
            Last Transcript:
          </div>
          <div style={{ fontSize: '0.7rem' }}>
            "{transcript.substring(0, 100)}{transcript.length > 100 ? '...' : ''}"
          </div>
        </div>
      )}

      {/* Error Display */}
      {error && (
        <div style={{
          marginBottom: '0.5rem',
          padding: '0.3rem',
          backgroundColor: 'rgba(255,0,0,0.1)',
          borderRadius: '0.25rem',
          border: '1px solid #ff6666'
        }}>
          <div style={{ color: '#ff6666', fontSize: '0.7rem', marginBottom: '0.2rem' }}>
            Current Error:
          </div>
          <div style={{ fontSize: '0.7rem', color: '#ff9999' }}>
            {error}
          </div>
        </div>
      )}

      {/* Controls */}
      <div style={{
        display: 'flex',
        gap: '0.5rem',
        flexWrap: 'wrap',
        marginTop: '0.5rem'
      }}>
        <button
          onClick={onForceStart}
          disabled={systemStatus.recognitionActive}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: systemStatus.recognitionActive ? '#666' : '#44ff44',
            color: '#000',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: systemStatus.recognitionActive ? 'not-allowed' : 'pointer',
            opacity: systemStatus.recognitionActive ? 0.5 : 1
          }}
        >
          🎤 Force Start
        </button>

        <button
          onClick={onForceStop}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: '#ff4444',
            color: '#fff',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: 'pointer'
          }}
        >
          🛑 Force Stop
        </button>

        <button
          onClick={onEmergencyReset}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: '#ff8800',
            color: '#fff',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: 'pointer',
            fontWeight: 'bold'
          }}
        >
          🚨 EMERGENCY RESET
        </button>
        
        <button
          onClick={() => {
            const controller = (window as any).voiceController
            if (controller) {
              // Manually activate conversation for testing
              controller.conversationActive = true
              controller.onConversationChange?.(true)
              console.log('🔧 Diagnostic: Manually activated conversation')
            }
          }}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: conversationActive ? '#666' : '#44ff44',
            color: conversationActive ? '#ccc' : '#000',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: 'pointer'
          }}
        >
          💬 Force Conversation
        </button>
        
        <button
          onClick={() => {
            const commandProcessor = (window as any).commandProcessor
            if (commandProcessor) {
              console.log('🔧 Diagnostic: Manually starting conversation via CommandProcessor')
              commandProcessor.setConversationActive(true)
            } else {
              console.error('🔧 Diagnostic: No CommandProcessor found')
            }
          }}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: '#00aa44',
            color: '#fff',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: 'pointer'
          }}
        >
          🎯 Force via Processor
        </button>
        
        <button
          onClick={() => {
            const socket = (window as any).socket || (document as any).socket
            if (socket) {
              console.log('🔧 Diagnostic: Sending test command to backend')
              socket.emit('voice_command', {
                command: 'test command from diagnostic panel',
                timestamp: Date.now(),
                conversationActive: true,
                isExplicitTrigger: false
              })
            } else {
              console.error('🔧 Diagnostic: No socket found')
            }
          }}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: '#0088ff',
            color: '#fff',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: 'pointer'
          }}
        >
          🧪 Test Backend
        </button>
        
        <button
          onClick={() => {
            const commandProcessor = (window as any).commandProcessor
            if (commandProcessor) {
              console.log('🔧 Diagnostic: Queuing test command via CommandProcessor')
              commandProcessor.queueCommand('test command from diagnostic', Date.now())
            } else {
              console.error('🔧 Diagnostic: No CommandProcessor found')
            }
          }}
          style={{
            fontSize: '0.6rem',
            padding: '0.3rem 0.6rem',
            backgroundColor: '#aa4400',
            color: '#fff',
            border: 'none',
            borderRadius: '0.25rem',
            cursor: 'pointer'
          }}
        >
          🔄 Queue Test
        </button>
      </div>

      {/* Detection Warnings */}
      {(isStuck || hasErrors || isInErrorState) && (
        <div style={{
          marginTop: '0.5rem',
          padding: '0.3rem',
          backgroundColor: 'rgba(255,165,0,0.2)',
          borderRadius: '0.25rem',
          border: '1px solid #ffaa00'
        }}>
          <div style={{ color: '#ffaa00', fontSize: '0.7rem', fontWeight: 'bold' }}>
            ⚠️ Issues Detected:
          </div>
          <div style={{ fontSize: '0.65rem', marginTop: '0.2rem' }}>
            {isStuck && <div>• System appears stuck (&gt;10s without transition)</div>}
            {hasErrors && <div>• {systemStatus.errorCount} error(s) occurred</div>}
            {isInErrorState && <div>• System is in ERROR state</div>}
          </div>
        </div>
      )}

      {/* Performance Metrics */}
      <div style={{
        marginTop: '0.5rem',
        fontSize: '0.65rem',
        opacity: 0.7,
        borderTop: '1px solid #444',
        paddingTop: '0.3rem'
      }}>
        <div>Controller Instance: {(window as any).voiceController ? '✅' : '❌'}</div>
        <div>WebAudio Support: {window.AudioContext ? '✅' : '❌'}</div>
        <div>Speech Recognition: {(window as any).webkitSpeechRecognition ? '✅' : '❌'}</div>
        <div>Page Focus: {document.hasFocus() ? '✅' : '❌'}</div>
      </div>

      {/* Instructions */}
      <div style={{
        marginTop: '0.5rem',
        fontSize: '0.6rem',
        opacity: 0.6,
        fontStyle: 'italic'
      }}>
        Click anywhere to update • This panel auto-refreshes every 100ms
      </div>

      {/* System Health Section */}
      <div style={{ marginBottom: '15px', padding: '10px', backgroundColor: '#1a1a1a', borderRadius: '6px' }}>
        <h4 style={{ margin: '0 0 8px 0', color: '#20808D', fontSize: '14px' }}>System Health</h4>
        <div style={{ fontSize: '12px', lineHeight: '1.4' }}>
          <div style={{ color: systemStatus.recognitionActive ? '#4ade80' : '#ef4444' }}>
            Recognition: {systemStatus.recognitionActive ? '✅' : '❌'}
          </div>
          <div style={{ color: systemStatus.ttsGenerating ? '#f59e0b' : systemStatus.ttsPlaying ? '#22c55e' : '#64748b' }}>
            TTS Gen: {systemStatus.ttsGenerating ? '❌' : '✅'} | TTS Play: {systemStatus.ttsPlaying ? '❌' : '✅'}
          </div>
          <div style={{ color: conversationActive ? '#4ade80' : '#64748b' }}>
            Conversation: {conversationActive ? '✅ ACTIVE' : '❌ INACTIVE'}
          </div>
          <div style={{ color: '#64748b' }}>
            Last Transition: {Math.round(timeSinceLastTransition / 1000)}s ago
          </div>
        </div>
      </div>

      {/* Conversation Persistence Info */}
      <div style={{ marginBottom: '15px', padding: '10px', backgroundColor: '#1a1a1a', borderRadius: '6px' }}>
        <h4 style={{ margin: '0 0 8px 0', color: '#20808D', fontSize: '14px' }}>Conversation Persistence</h4>
        <div style={{ fontSize: '12px', lineHeight: '1.4' }}>
          <div style={{ color: conversationActive ? '#4ade80' : '#ef4444' }}>
            State: {conversationActive ? 'MAINTAINED' : 'RESET'}
          </div>
          <div style={{ color: '#64748b' }}>
            Follow-ups: {conversationActive ? 'ENABLED' : 'DISABLED'}
          </div>
          <div style={{ color: '#64748b' }}>
            Error Count: {systemStatus.errorCount}
          </div>
        </div>
      </div>
    </div>
  )
}

export default DiagnosticOverlay 