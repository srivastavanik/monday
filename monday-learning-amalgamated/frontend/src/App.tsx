import React, { useEffect, useState, useCallback } from 'react'
import { Canvas } from '@react-three/fiber'
import { XR, Controllers, Hands } from '@react-three/xr'
import MondayScene from './components/MondayScene'
import VoiceInterface from './components/VoiceInterface'
import SpatialOrchestrator from './components/SpatialOrchestrator'
import LoadingOverlay from './components/LoadingOverlay'
import ErrorBoundary from './components/ErrorBoundary'
import DiagnosticOverlay from './components/DiagnosticOverlay'
import StaticInfoPanel from './components/StaticInfoPanel'
import { useMondayStore } from './store/mondayStore'
import { useVoiceSystem } from './hooks/useVoiceSystem'
import { useWebSocketConnection } from './hooks/useWebSocket'
import { usePerformanceMonitor } from './hooks/usePerformanceMonitor'
import { SystemState } from './controllers/VoiceSystemController'

const App: React.FC = () => {
  const [isInitialized, setIsInitialized] = useState(false)
  const [isVRSupported, setIsVRSupported] = useState(false)
  const [lastProcessedTranscript, setLastProcessedTranscript] = useState('')
  const [audioInitialized, setAudioInitialized] = useState(false)
  const [currentModel, setCurrentModel] = useState<string>('')
  const [isProcessingReasoning, setIsProcessingReasoning] = useState(false)
  const [staticPanelData, setStaticPanelData] = useState<{
    isVisible: boolean
    title: string
    content: string
    citations: any[]
    model: string
  }>({
    isVisible: false,
    title: '',
    content: '',
    citations: [],
    model: ''
  })

  const { 
    isConnected, 
    sessionState, 
    initializeSession,
    setPerformanceMetrics,
    addPanel,
    setActivePanel,
    setConversationActive,
    updatePanel,
    panels
  } = useMondayStore()

  // Initialize WebXR and check support
  useEffect(() => {
    const checkVRSupport = async () => {
      if (navigator.xr) {
        try {
          const supported = await navigator.xr.isSessionSupported('immersive-vr')
          setIsVRSupported(supported)
          console.log('VR Support:', supported ? 'Available' : 'Not available')
        } catch (error) {
          console.warn('VR support check failed:', error)
          setIsVRSupported(false)
        }
      } else {
        console.warn('WebXR not available')
        setIsVRSupported(false)
      }
      setIsInitialized(true)
    }

    checkVRSupport()
  }, [])

  // Initialize WebSocket connection
  const { socket, isConnected: socketConnected } = useWebSocketConnection()

  // Initialize the unified voice system
  const {
    systemState,
    isListening,
    isSpeaking,
    isPlaying,
    transcript,
    conversationActive,
    error: voiceError,
    systemStatus,
    handleCommand,
    handleTTSResponse,
    emergencyReset,
    forceStartListening,
    forceStop,
    interruptTTS,
    lockMicrophoneForProcess,
    unlockMicrophoneAfterProcess
  } = useVoiceSystem({
    onError: (error) => {
      console.error('🌟 App: Voice system error:', error)
    }
  })

  // Initialize audio on first user interaction
  const handleUserInteraction = useCallback(async (): Promise<boolean> => {
    if (!audioInitialized) {
      try {
        // Create a temporary audio context to test
        const tempContext = new (window.AudioContext || (window as any).webkitAudioContext)()
        if (tempContext.state === 'suspended') {
          await tempContext.resume()
        }
        await tempContext.close()
        
        setAudioInitialized(true)
        console.log('🌟 App: Audio initialized after user interaction')
        return true
      } catch (error) {
        console.warn('🌟 App: Audio initialization failed:', error)
        return false
      }
    }
    return true
  }, [audioInitialized])

  // Handle interrupting TTS
  const handleInterruptTTS = useCallback(async () => {
    if (systemState === SystemState.PLAYING_TTS) {
      console.log('🌟 App: Interrupting TTS and starting listening')
      await interruptTTS() // Use the new interrupt method
    } else {
      await forceStartListening()
    }
  }, [systemState, interruptTTS, forceStartListening])

  // Handle WebSocket responses from backend
  useEffect(() => {
    if (!socket) return

    // Make socket globally available for ProgressivePanel
    ;(window as any).socket = socket
    console.log('🌟 App: Socket made globally available for progressive panels')

    const handleVoiceResponse = async (response: any) => {
      console.log('🌟 App: Received voice response from backend:', response)
      
      // Update current model indicator
      if (response.data?.model) {
        setCurrentModel(response.data.model)
      }
      
      // Update static panel with full response content
      if (response.data?.panels && Array.isArray(response.data.panels) && response.data.panels.length > 0) {
        const mainPanel = response.data.panels[0] // Use the first panel as the main content
        
        // For reasoning/research responses, show the full content, not the TTS message
        let displayContent = mainPanel.content || 'No content available'
        let displayTitle = mainPanel.title || 'Monday Response'
        
        // If there's fullContent available (reasoning/research), use that instead
        if (mainPanel.fullContent && mainPanel.fullContent !== mainPanel.content) {
          displayContent = mainPanel.fullContent
          
          // Clean up the reasoning content for display
          if (displayContent.startsWith('<think>')) {
            // Extract the reasoning content from the <think> tags
            const thinkMatch = displayContent.match(/<think>([\s\S]*?)<\/think>/);
            if (thinkMatch) {
              displayContent = thinkMatch[1].trim()
            } else {
              // Remove <think> tags if no closing tag
              displayContent = displayContent.replace(/<think>\s*/g, '').trim()
            }
          }
          
          // Update title for reasoning/research
          if (response.data.mode === 'reasoning') {
            displayTitle = 'Monday\'s Reasoning Process'
          } else if (response.data.mode === 'research') {
            displayTitle = 'Monday\'s Research Analysis'
          }
        }
        
        console.log('🌟 App: Setting static panel data:', {
          title: displayTitle,
          contentLength: displayContent.length,
          contentPreview: displayContent.substring(0, 100),
          mode: response.data.mode,
          hasFullContent: !!mainPanel.fullContent
        })
        
        setStaticPanelData({
          isVisible: true,
          title: displayTitle,
          content: displayContent,
          citations: response.data.citations || mainPanel.citations || [],
          model: response.data.model || 'sonar'
        })
      } else {
        // Fallback: create static panel from response data directly
        let displayContent = response.message || 'No content available'
        let displayTitle = 'Monday Response'
        
        // Check if this is a reasoning/research response with metadata
        if (response.data?.metadata?.isThinking && response.data?.metadata?.fullContent) {
          displayContent = response.data.metadata.fullContent
          displayTitle = 'Monday\'s Reasoning Process'
          
          // Clean up reasoning content
          if (displayContent.startsWith('<think>')) {
            const thinkMatch = displayContent.match(/<think>([\s\S]*?)<\/think>/);
            if (thinkMatch) {
              displayContent = thinkMatch[1].trim()
            } else {
              displayContent = displayContent.replace(/<think>\s*/g, '').trim()
            }
          }
        } else if (response.data?.metadata?.isResearching && response.data?.metadata?.fullContent) {
          displayContent = response.data.metadata.fullContent
          displayTitle = 'Monday\'s Research Analysis'
        }
        
        console.log('🌟 App: Setting fallback static panel data:', {
          title: displayTitle,
          contentLength: displayContent.length,
          contentPreview: displayContent.substring(0, 100),
          mode: response.data?.mode,
          isThinking: response.data?.metadata?.isThinking,
          isResearching: response.data?.metadata?.isResearching
        })
        
        setStaticPanelData({
          isVisible: true,
          title: displayTitle,
          content: displayContent,
          citations: response.data?.citations || [],
          model: response.data?.model || 'sonar'
        })
      }
      
      // Add panels to the store for 3D visualization
      if (response.data?.panels && Array.isArray(response.data.panels)) {
        response.data.panels.forEach((panelData: any) => {
          console.log('🌟 App: Adding panel to store:', panelData)
          addPanel(panelData)
        })
        
        // Set the first panel as active if there are panels
        if (response.data.panels.length > 0) {
          setActivePanel(response.data.panels[0].id)
        }
      }
      
      // Handle TTS response through the unified system
      if (response.message) {
        console.log('🌟 App: Triggering TTS through voice system controller')
        try {
          await handleTTSResponse(response.message)
          console.log('🌟 App: TTS handling completed successfully')
        } catch (error) {
          console.error('🌟 App: TTS handling failed:', error)
        }
      }
    }

    const handleVoiceError = (errorData: any) => {
      console.error('🌟 App: Backend voice error:', errorData)
    }

    const handleReasoningProgress = async (data: any) => {
      console.log('🔄 App: Reasoning progress update:', data)
      
      // Set processing flag when reasoning starts
      if (data.update?.progress === 10 && !isProcessingReasoning) {
        setIsProcessingReasoning(true)
        console.log('🧠 App: Started reasoning process')
      }
      
      // Update static panel with real-time reasoning content (but don't trigger TTS during progress)
      if (data.update?.reasoning && data.update?.type !== 'complete') {
        let cleanContent = data.update.reasoning
        
        // Clean up the reasoning content for display
        if (cleanContent.startsWith('<think>')) {
          const thinkMatch = cleanContent.match(/<think>([\s\S]*?)<\/think>/);
          if (thinkMatch) {
            cleanContent = thinkMatch[1].trim()
          } else {
            // Remove opening <think> tag if no closing tag yet (streaming in progress)
            cleanContent = cleanContent.replace(/^<think>\s*/g, '').trim()
          }
        }
        
        // PARAGRAPH-LEVEL ROTATION: Show the latest complete paragraph
        let displayContent = cleanContent
        if (cleanContent.length > 0) {
          // Split into paragraphs (double newlines or single newlines with substantial content)
          const paragraphs = cleanContent.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 50)
          
          if (paragraphs.length > 0) {
            // Show the latest paragraph, or last 2 paragraphs if they're short
            if (paragraphs.length === 1) {
              displayContent = paragraphs[0].trim()
            } else if (paragraphs.length >= 2) {
              const lastTwo = paragraphs.slice(-2).join('\n\n').trim()
              // If the combined length is reasonable, show both, otherwise just the latest
              displayContent = lastTwo.length > 800 ? paragraphs[paragraphs.length - 1].trim() : lastTwo
            }
            
            // Add indicator if there are more paragraphs
            if (paragraphs.length > 2) {
              displayContent = `[Continuing analysis...]\n\n${displayContent}`
            }
          }
        }
        
        // Ensure we have meaningful content to display
        if (displayContent.length > 10) {
          // Update static panel with paragraph-level reasoning content
          setStaticPanelData({
            isVisible: true,
            title: `Monday's Reasoning Process (${data.update.progress || 0}% complete)`,
            content: displayContent,
            citations: [],
            model: data.model || 'sonar-reasoning-pro'
          })
          
          console.log('🔄 App: Updated static panel with paragraph-level reasoning:', {
            progress: data.update.progress,
            originalLength: cleanContent.length,
            paragraphCount: cleanContent.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 50).length,
            displayLength: displayContent.length,
            isComplete: data.update.type === 'complete'
          })
        }
      }
      
      // Update the reasoning panel in the store if it exists
      if (data.update?.reasoning && panels.length >= 2) {
        const reasoningPanel = panels.find(p => p.type === 'reasoning')
        if (reasoningPanel) {
          let cleanContent = data.update.reasoning
          if (cleanContent.startsWith('<think>')) {
            const thinkMatch = cleanContent.match(/<think>([\s\S]*?)<\/think>/);
            if (thinkMatch) {
              cleanContent = thinkMatch[1].trim()
            } else {
              cleanContent = cleanContent.replace(/^<think>\s*/g, '').trim()
            }
          }
          
          updatePanel(reasoningPanel.id, {
            content: cleanContent,
            fullContent: data.update.reasoning,
            title: `Reasoning Process (${data.update.progress || 0}% complete)`
          })
        }
      }
      
      // When reasoning is complete, show FULL FINAL RESPONSE and provide substantial TTS (ONLY ONCE)
      if (data.update?.type === 'complete' && data.update?.progress === 100 && isProcessingReasoning) {
        console.log('🔓 App: Reasoning completed, displaying full final response')
        setIsProcessingReasoning(false) // Reset processing flag
        
        // Get the complete final reasoning content
        let finalContent = data.update.reasoning || ''
        if (finalContent.startsWith('<think>')) {
          const thinkMatch = finalContent.match(/<think>([\s\S]*?)<\/think>/);
          if (thinkMatch) {
            finalContent = thinkMatch[1].trim()
          } else {
            finalContent = finalContent.replace(/^<think>\s*/g, '').trim()
          }
        }
        
        // Display the COMPLETE final response in the static panel
        setStaticPanelData({
          isVisible: true,
          title: 'Monday\'s Complete Reasoning Analysis ✅',
          content: finalContent, // FULL CONTENT, not a summary
          citations: [],
          model: data.model || 'sonar-reasoning-pro'
        })
        
        // Create substantial TTS content that explains the key insights
        let substantialTTS = ''
        if (finalContent.length > 100) {
          // Extract the main points for a substantial spoken explanation
          const paragraphs = finalContent.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 30)
          
          if (paragraphs.length >= 2) {
            // Create a substantial explanation from the first few paragraphs
            const keyInsights = paragraphs.slice(0, 2).join(' ').trim()
            const firstSentences = keyInsights.split(/[.!?]+/).slice(0, 4).join('. ').trim()
            substantialTTS = `I've completed my reasoning analysis. Here's what I found: ${firstSentences}. The complete analysis is now displayed for you to explore in detail.`
          } else if (paragraphs.length === 1) {
            // Use the single paragraph but make it more conversational
            const sentences = paragraphs[0].split(/[.!?]+/).filter((s: string) => s.trim().length > 10).slice(0, 3)
            substantialTTS = `I've finished thinking through this step by step. ${sentences.join('. ')}. You can see my complete reasoning process in the panel.`
          } else {
            // Fallback for edge cases
            const sentences = finalContent.split(/[.!?]+/).filter((s: string) => s.trim().length > 15).slice(0, 3)
            substantialTTS = `I've completed my analysis. ${sentences.join('. ')}. The full reasoning is available for you to review.`
          }
        } else {
          substantialTTS = `I've completed my reasoning analysis. The full thought process is now available for you to explore in the panel.`
        }
        
        console.log('🎯 App: Displaying complete final response with substantial TTS:', {
          finalContentLength: finalContent.length,
          ttsLength: substantialTTS.length,
          ttsPreview: substantialTTS.substring(0, 100)
        })
        
        // Trigger substantial TTS (ONLY ONCE)
        try {
          await handleTTSResponse(substantialTTS)
        } catch (error) {
          console.error('🔊 App: Failed to play substantial TTS:', error)
        }
        
        // Unlock microphone after a delay to allow TTS to complete
        setTimeout(async () => {
          try {
            await unlockMicrophoneAfterProcess()
          } catch (error) {
            console.error('🔓 App: Failed to unlock microphone after reasoning:', error)
          }
        }, 4000) // Longer delay for substantial TTS
      }
    }

    const handleResearchProgress = async (data: any) => {
      console.log('🔍 App: Research progress update:', data)
      
      // Update static panel with real-time research content
      if (data.update?.reasoning) {
        // PARAGRAPH-LEVEL ROTATION: Show the latest complete paragraph
        let displayContent = data.update.reasoning
        if (data.update.reasoning.length > 0) {
          // Split into paragraphs (double newlines or single newlines with substantial content)
          const paragraphs = data.update.reasoning.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 50)
          
          if (paragraphs.length > 0) {
            // Show the latest paragraph, or last 2 paragraphs if they're short
            if (paragraphs.length === 1) {
              displayContent = paragraphs[0].trim()
            } else if (paragraphs.length >= 2) {
              const lastTwo = paragraphs.slice(-2).join('\n\n').trim()
              // If the combined length is reasonable, show both, otherwise just the latest
              displayContent = lastTwo.length > 800 ? paragraphs[paragraphs.length - 1].trim() : lastTwo
            }
            
            // Add indicator if there are more paragraphs
            if (paragraphs.length > 2) {
              displayContent = `[Continuing research...]\n\n${displayContent}`
            }
          }
        }
        
        // Update static panel with paragraph-level research content
        setStaticPanelData({
          isVisible: true,
          title: `Monday's Research Analysis (${data.update.progress || 0}% complete)`,
          content: displayContent,
          citations: data.update?.sources || [],
          model: data.model || 'sonar-deep-research'
        })
        
        console.log('🔍 App: Updated static panel with paragraph-level research:', {
          progress: data.update.progress,
          originalLength: data.update.reasoning.length,
          paragraphCount: data.update.reasoning.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 50).length,
          displayLength: displayContent.length,
          isComplete: data.update.type === 'complete'
        })
      }
      
      // Update the research panel in the store if it exists
      if (data.update?.reasoning && panels.length >= 2) {
        const researchPanel = panels.find(p => p.type === 'reasoning')
        if (researchPanel) {
          updatePanel(researchPanel.id, {
            content: data.update.reasoning,
            fullContent: data.update.reasoning,
            title: `Research Analysis (${data.update.progress || 0}% complete)`
          })
        }
      }
      
      // When research is complete, show FULL FINAL RESPONSE and provide substantial TTS
      if (data.update?.type === 'complete' && data.update?.progress === 100) {
        console.log('🔓 App: Research completed, displaying full final response')
        
        // Get the complete final research content
        const finalContent = data.update?.reasoning || ''
        
        // Display the COMPLETE final response in the static panel
        setStaticPanelData({
          isVisible: true,
          title: 'Monday's Complete Research Analysis ✅',
          content: finalContent, // FULL CONTENT, not a summary
          citations: data.update?.sources || [],
          model: data.model || 'sonar-deep-research'
        })
        
        // Create substantial TTS content that explains the key research findings
        let substantialTTS = ''
        if (finalContent.length > 100) {
          // Extract the main research findings for a substantial spoken explanation
          const paragraphs = finalContent.split(/\n\s*\n|\n(?=[A-Z])/g).filter((p: string) => p.trim().length > 30)
          
          if (paragraphs.length >= 2) {
            // Create a substantial explanation from the first few paragraphs
            const keyFindings = paragraphs.slice(0, 2).join(' ').trim()
            const firstSentences = keyFindings.split(/[.!?]+/).slice(0, 4).join('. ').trim()
            substantialTTS = `I've completed my comprehensive research analysis. Here are the key findings: ${firstSentences}. The complete research with sources is now displayed for you to explore in detail.`
          } else if (paragraphs.length === 1) {
            // Use the single paragraph but make it more conversational
            const sentences = paragraphs[0].split(/[.!?]+/).filter((s: string) => s.trim().length > 10).slice(0, 3)
            substantialTTS = `I've finished researching this topic thoroughly. ${sentences.join('. ')}. You can see my complete research analysis with sources in the panel.`
          } else {
            // Fallback for edge cases
            const sentences = finalContent.split(/[.!?]+/).filter((s: string) => s.trim().length > 15).slice(0, 3)
            substantialTTS = `I've completed my research analysis. ${sentences.join('. ')}. The full research findings with sources are available for you to review.`
          }
        } else {
          substantialTTS = `I've completed my research analysis. The full research findings with sources are now available for you to explore in the panel.`
        }
        
        console.log('🎯 App: Displaying complete research response with substantial TTS:', {
          finalContentLength: finalContent.length,
          ttsLength: substantialTTS.length,
          ttsPreview: substantialTTS.substring(0, 100)
        })
        
        // Trigger substantial TTS
        try {
          await handleTTSResponse(substantialTTS)
        } catch (error) {
          console.error('🔊 App: Failed to play substantial research TTS:', error)
        }
        
        // Unlock microphone after a delay to allow TTS to complete
        setTimeout(async () => {
          try {
            await unlockMicrophoneAfterProcess()
          } catch (error) {
            console.error('🔓 App: Failed to unlock microphone after research:', error)
          }
        }, 3000) // Longer delay for substantial TTS
      }
    }

    socket.on('voice_response', handleVoiceResponse)
    socket.on('voice_error', handleVoiceError)
    socket.on('reasoning_progress', handleReasoningProgress)
    socket.on('research_progress', handleResearchProgress)

    return () => {
      socket.off('voice_response', handleVoiceResponse)
      socket.off('voice_error', handleVoiceError)
      socket.off('reasoning_progress', handleReasoningProgress)
      socket.off('research_progress', handleResearchProgress)
    }
  }, [socket, handleTTSResponse, addPanel, setActivePanel, updatePanel, panels, unlockMicrophoneAfterProcess])

  // Performance monitoring
  const { fps, frameTime, memoryUsage } = usePerformanceMonitor()

  // Update performance metrics in store
  useEffect(() => {
    setPerformanceMetrics({
      fps,
      frameTime,
      memoryUsage,
      timestamp: Date.now()
    })
  }, [fps, frameTime, memoryUsage, setPerformanceMetrics])

  // Initialize Monday session
  useEffect(() => {
    if (isInitialized && socketConnected) {
      initializeSession({
        isVRSupported,
        userAgent: navigator.userAgent,
        timestamp: Date.now()
      })
    }
  }, [isInitialized, socketConnected, isVRSupported, initializeSession])

  // Process voice commands when transcript changes (avoids stale closure issues)
  useEffect(() => {
    if (!transcript || transcript === lastProcessedTranscript) {
      return
    }

    console.log('🌟 App: New transcript detected, processing command:', {
      transcript: transcript,
      lastProcessed: lastProcessedTranscript,
      socketConnected,
      conversationActive,
      hasSocket: !!socket,
      socketId: socket?.id || 'no-socket'
    })

    // Check prerequisites with current values (not stale closures)
    if (!socket || !socketConnected || !transcript.trim()) {
      console.log('🌟 App: Command processing blocked - missing prerequisites:', {
        hasSocket: !!socket,
        socketConnected,
        hasCommand: !!transcript.trim(),
        socketReadyState: socket?.connected
      })
      return
    }

    const normalizedCommand = transcript.toLowerCase().trim()
    
    // Check if this should trigger Monday (either explicit trigger or in conversation)
    const isExplicitTrigger = normalizedCommand.includes('hey monday')
    const shouldProcess = isExplicitTrigger || conversationActive
    
    // CRITICAL: Filter out Monday's own voice to prevent feedback loops
    const isMondayVoice = normalizedCommand.includes('think through') || 
                         normalizedCommand.includes('machine learning algorithms step by step') ||
                         normalizedCommand.includes('break down think about') ||
                         normalizedCommand.includes('analyze the different aspects') ||
                         normalizedCommand.includes('completed my reasoning analysis') ||
                         normalizedCommand.includes('here\'s what i found')
    
    if (isMondayVoice) {
      console.log('🚫 App: Filtering out Monday\'s own voice to prevent feedback loop:', {
        command: normalizedCommand.substring(0, 50)
      })
      return
    }
    
    console.log('🌟 App: Command filtering:', {
      command: normalizedCommand.substring(0, 50),
      isExplicitTrigger,
      conversationActive,
      shouldProcess,
      isMondayVoice
    })

    if (shouldProcess) {
      // CRITICAL: Lock microphone IMMEDIATELY for reasoning/research commands
      const isReasoningCommand = normalizedCommand.includes('think') || normalizedCommand.includes('reason')
      const isResearchCommand = normalizedCommand.includes('research') || normalizedCommand.includes('investigate')
      
      if (isReasoningCommand) {
        console.log('🔇 App: IMMEDIATELY locking microphone for reasoning command')
        lockMicrophoneForProcess('reasoning').catch(error => {
          console.error('🔇 App: Failed to lock microphone for reasoning:', error)
        })
      } else if (isResearchCommand) {
        console.log('🔇 App: IMMEDIATELY locking microphone for research command')
        lockMicrophoneForProcess('research').catch(error => {
          console.error('🔇 App: Failed to lock microphone for research:', error)
        })
      }
      
      // If it's an explicit trigger, set conversation active
      if (isExplicitTrigger) {
        setConversationActive(true)
      }

      console.log('🌟 App: 🎯 Processing voice command:', {
        command: transcript,
        isExplicitTrigger,
        conversationActive,
        commandLength: transcript.length,
        isReasoningCommand,
        isResearchCommand
      })
      
      console.log('🌟 App: 📤 Sending command to backend via WebSocket')
      
      try {
        socket.emit('voice_command', {
          command: transcript,
          timestamp: Date.now(),
          conversationActive: true, // Always send true if we're processing
          isExplicitTrigger: isExplicitTrigger
        })
        
        setLastProcessedTranscript(transcript)
        console.log('🌟 App: ✅ Command sent to backend successfully')
      } catch (error) {
        console.error('🌟 App: ❌ Failed to send command to backend:', error)
      }
    } else {
      console.log('🌟 App: ❌ Command ignored - no trigger and not in conversation:', {
        command: normalizedCommand.substring(0, 50),
        needsTrigger: !isExplicitTrigger && !conversationActive
      })
    }
  }, [transcript, lastProcessedTranscript, socket, socketConnected, conversationActive, lockMicrophoneForProcess])

  if (!isInitialized) {
    return <LoadingOverlay message="Initializing Monday..." />
  }

  // Check for critical voice errors that require full-screen handling
  if (voiceError && voiceError.includes('Microphone permission lost')) {
    return (
      <div style={{ 
        padding: '2rem', 
        textAlign: 'center', 
        color: 'var(--paper-white)',
        backgroundColor: 'var(--offblack)',
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center'
      }}>
        <h2>Microphone Access Required</h2>
        <p>Monday needs microphone access to function properly.</p>
        <p>Please refresh the page and allow microphone access when prompted.</p>
        <button 
          className="btn" 
          onClick={() => window.location.reload()}
          style={{ margin: '1rem auto', display: 'block' }}
        >
          Refresh Page
        </button>
      </div>
    )
  }

  return (
    <ErrorBoundary>
      <div style={{ 
        width: '100vw', 
        height: '100vh', 
        backgroundColor: 'var(--offblack)' 
      }}>
        <Canvas
          camera={{ 
            position: [0, 1.6, 0],
            fov: 75 
          }}
          gl={{ 
            antialias: true,
            alpha: false,
            powerPreference: 'high-performance' 
          }}
          frameloop="always"
        >
          <XR referenceSpace="local-floor">
            <ambientLight intensity={0.3} color="#20808D" />
            <directionalLight 
              position={[5, 5, 5]} 
              intensity={0.5} 
              color="#FBFAF4" 
            />
            
            <Controllers />
            <Hands />
            <MondayScene />
            <SpatialOrchestrator />
          </XR>
        </Canvas>

        <VoiceInterface 
          isListening={isListening}
          transcript={transcript}
          onStartListening={handleInterruptTTS}
          onStopListening={forceStop}
          conversationActive={conversationActive}
          onUserInteraction={handleUserInteraction}
          isRestartingVoice={systemState === SystemState.RESETTING}
          isSpeaking={isSpeaking}
          isPlaying={isPlaying}
        />

        {/* Model Indicator */}
        {currentModel && (
          <div style={{
            position: 'fixed',
            top: '1rem',
            right: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'rgba(32, 128, 141, 0.9)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            border: '1px solid var(--true-turquoise)'
          }}>
            🤖 {currentModel.replace('sonar-', '').replace('-', ' ').toUpperCase()}
          </div>
        )}

        {/* Diagnostic Overlay - Only on localhost */}
        {window.location.hostname === 'localhost' && (
          <DiagnosticOverlay
            systemState={systemState}
            systemStatus={systemStatus}
            transcript={transcript}
            conversationActive={conversationActive}
            error={voiceError}
            onEmergencyReset={emergencyReset}
            onForceStart={forceStartListening}
            onForceStop={forceStop}
          />
        )}

        {/* Connection Status */}
        <div style={{
          position: 'fixed',
          bottom: '1rem',
          right: '1rem',
          padding: '0.5rem 1rem',
          backgroundColor: socketConnected ? 'var(--true-turquoise)' : '#dc3545',
          color: 'var(--paper-white)',
          borderRadius: '0.25rem',
          fontSize: '0.875rem',
          zIndex: 1000
        }}>
          {socketConnected ? 'Connected' : 'Disconnected'}
        </div>

        {/* Debug: Panel Count */}
        {window.location.hostname === 'localhost' && (
          <div style={{
            position: 'fixed',
            bottom: '1rem',
            right: '8rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'rgba(32, 128, 141, 0.8)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000
          }}>
            Panels: {panels.length}
          </div>
        )}

        {/* Debug: Panel Details */}
        {window.location.hostname === 'localhost' && panels.length > 0 && (
          <div style={{
            position: 'fixed',
            top: '1rem',
            left: '1rem',
            padding: '1rem',
            backgroundColor: 'rgba(9, 23, 23, 0.9)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.75rem',
            zIndex: 1000,
            maxWidth: '300px',
            maxHeight: '200px',
            overflow: 'auto'
          }}>
            <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
              Debug: Panels in Store ({panels.length})
            </div>
            {panels.map((panel, index) => (
              <div key={panel.id} style={{ marginBottom: '0.5rem', padding: '0.25rem', backgroundColor: 'rgba(32, 128, 141, 0.2)' }}>
                <div>#{index + 1}: {panel.title}</div>
                <div>Type: {panel.type}</div>
                <div>ID: {panel.id}</div>
                <div>Active: {panel.isActive ? 'Yes' : 'No'}</div>
                <div>Position: [{panel.position.join(', ')}]</div>
              </div>
            ))}
          </div>
        )}

        {/* Conversation Status */}
        {conversationActive && (
          <div style={{
            position: 'fixed',
            bottom: '1rem',
            left: '1rem',
            padding: '0.5rem 1rem',
            backgroundColor: 'var(--true-turquoise)',
            color: 'var(--paper-white)',
            borderRadius: '0.25rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            display: 'flex',
            alignItems: 'center',
            gap: '0.5rem'
          }}>
            <div style={{
              width: '8px',
              height: '8px',
              backgroundColor: 'var(--paper-white)',
              borderRadius: '50%',
              animation: 'pulse 2s ease-in-out infinite'
            }} />
            Conversation Active
          </div>
        )}

        {/* System Status */}
        <div style={{
          position: 'fixed',
          bottom: '3.5rem',
          right: '1rem',
          padding: '0.5rem 1rem',
          backgroundColor: systemState === SystemState.ERROR ? '#dc3545' : 
                          systemState === SystemState.ACTIVE_LISTENING ? 'var(--true-turquoise)' :
                          systemState === SystemState.PLAYING_TTS ? '#ff6600' :
                          'rgba(9, 23, 23, 0.8)',
          color: 'var(--paper-white)',
          borderRadius: '0.25rem',
          fontSize: '0.875rem',
          zIndex: 1000,
          display: 'flex',
          alignItems: 'center',
          gap: '0.5rem'
        }}>
          {systemState === SystemState.PLAYING_TTS && (
            <div style={{
              width: '8px',
              height: '8px',
              backgroundColor: 'var(--paper-white)',
              borderRadius: '50%',
              animation: 'pulse 1s ease-in-out infinite'
            }} />
          )}
          {systemState === SystemState.ACTIVE_LISTENING ? '🎤 Listening' :
           systemState === SystemState.PLAYING_TTS ? '🔊 Speaking (Mic Muted)' :
           systemState === SystemState.PROCESSING_COMMAND ? '⚙️ Processing' :
           systemState === SystemState.RESETTING ? '🔄 Resetting' :
           systemState === SystemState.ERROR ? '❌ Error' :
           systemState === SystemState.WAITING_FOR_ACTIVATION ? '⏳ Waiting' :
           '💤 Idle'}
        </div>

        {/* Audio Initialization Prompt */}
        {!audioInitialized && (
          <div 
            onClick={handleUserInteraction}
            style={{
              position: 'fixed',
              bottom: '6rem',
              right: '1rem',
              padding: '1rem',
              backgroundColor: '#f39c12',
              color: 'var(--paper-white)',
              borderRadius: '0.5rem',
              fontSize: '0.875rem',
              zIndex: 1000,
              cursor: 'pointer',
              border: '2px solid #e67e22'
            }}
          >
            🔇 Click to initialize audio for voice responses
          </div>
        )}

        {/* Error Display */}
        {voiceError && !voiceError.includes('Microphone permission lost') && (
          <div style={{
            position: 'fixed',
            bottom: '8rem',
            right: '1rem',
            padding: '1rem',
            backgroundColor: '#dc3545',
            color: 'var(--paper-white)',
            borderRadius: '0.5rem',
            fontSize: '0.875rem',
            zIndex: 1000,
            maxWidth: '400px',
            cursor: 'pointer'
          }}
          onClick={emergencyReset}
          >
            <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
               Voice System Error
            </div>
            <div style={{ fontSize: '0.75rem', marginBottom: '0.5rem' }}>
              {voiceError}
            </div>
            <div style={{ fontSize: '0.7rem', opacity: 0.8 }}>
              Click to reset the voice system.
            </div>
          </div>
        )}

        {/* Welcome Instructions */}
        {!sessionState.hasCompletedIntro && !conversationActive && systemState === SystemState.WAITING_FOR_ACTIVATION && (
          <div style={{
            position: 'fixed',
            top: '50%',
            left: '50%',
            transform: 'translate(-50%, -50%)',
            padding: '2rem',
            backgroundColor: 'var(--paper-white)',
            color: 'var(--offblack)',
            borderRadius: '0.5rem',
            border: '2px solid var(--true-turquoise)',
            textAlign: 'center',
            maxWidth: '500px',
            zIndex: 2000
          }}>
            <h3 style={{ color: 'var(--true-turquoise)', marginBottom: '1rem' }}>
              Welcome to Monday
            </h3>
            <p><strong>Say "Hey Monday"</strong> to start your learning journey.</p>
            <p style={{ fontSize: '0.875rem', margin: '1rem 0' }}>
              After that, you can ask questions without saying "Hey Monday" each time:
            </p>
            <ul style={{ fontSize: '0.875rem', textAlign: 'left', color: 'var(--true-turquoise)' }}>
              <li>"Hey Monday, tell me about quantum physics"</li>
              <li>"How does it work?" (no Hey Monday needed)</li>
              <li>"Think about machine learning algorithms"</li>
              <li>"Research the latest developments"</li>
            </ul>
            <div style={{
              fontSize: '0.75rem',
              marginTop: '1rem',
              padding: '0.5rem',
              backgroundColor: 'rgba(32, 128, 141, 0.1)',
              borderRadius: '0.25rem'
            }}>
              <strong>System Status:</strong> {systemState} | 
              <strong> Voice:</strong> {isListening ? ' 🎤 Ready' : ' 🔇 Not Active'}
            </div>
          </div>
        )}

        {/* Global Watchdog Alert */}
        {window.location.hostname === 'localhost' && systemStatus.errorCount > 3 && (
          <div 
            onClick={emergencyReset}
            style={{
              position: 'fixed',
              top: '50%',
              right: '1rem',
              transform: 'translateY(-50%)',
              padding: '1rem',
              backgroundColor: '#ff4500',
              color: 'var(--paper-white)',
              borderRadius: '0.5rem',
              fontSize: '0.875rem',
              zIndex: 2001,
              cursor: 'pointer',
              border: '3px solid #ff6500',
              maxWidth: '300px',
              fontWeight: 'bold'
            }}
          >
            🚨 SYSTEM UNSTABLE
            <div style={{ fontSize: '0.7rem', marginTop: '0.5rem', fontWeight: 'normal' }}>
              {systemStatus.errorCount} errors detected. Click for emergency reset.
            </div>
          </div>
        )}

        {/* Microphone Muted Indicator */}
        {systemState === SystemState.PLAYING_TTS && (
          <div style={{
            position: 'fixed',
            top: '50%',
            left: '20px',
            transform: 'translateY(-50%)',
            padding: '1rem',
            backgroundColor: '#ff6600',
            color: 'var(--paper-white)',
            borderRadius: '0.5rem',
            fontSize: '0.875rem',
            zIndex: 1001,
            display: 'flex',
            alignItems: 'center',
            gap: '0.5rem',
            border: '2px solid #ff8800',
            boxShadow: '0 4px 12px rgba(255, 102, 0, 0.3)'
          }}>
            <div style={{
              width: '12px',
              height: '12px',
              backgroundColor: 'var(--paper-white)',
              borderRadius: '50%',
              animation: 'pulse 1s ease-in-out infinite'
            }} />
            <div>
              <div style={{ fontWeight: 'bold' }}>🔇 Microphone Muted</div>
              <div style={{ fontSize: '0.75rem', opacity: 0.9 }}>
                Preventing feedback during speech
              </div>
            </div>
          </div>
        )}

        {/* Static Information Panel */}
        <StaticInfoPanel
          id="main-response"
          title={staticPanelData.title}
          content={staticPanelData.content}
          citations={staticPanelData.citations}
          model={staticPanelData.model}
          isVisible={staticPanelData.isVisible}
          onClose={() => setStaticPanelData(prev => ({ ...prev, isVisible: false }))}
        />
      </div>
    </ErrorBoundary>
  )
}

export default App 