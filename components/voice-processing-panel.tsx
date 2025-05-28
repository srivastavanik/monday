"use client"

import { useState, useEffect } from "react"
import { Mic, Volume2, Brain, Zap, Activity, Sparkles, MessageSquare, Loader2 } from "lucide-react"
import { PerplexityLogo } from "./perplexity-logo"

interface VoiceProcessingPanelProps {
  currentResponse: string
  isThinking: boolean
  isTyping: boolean
  mode: "basic" | "reasoning" | "deep-research"
  userTranscript?: string
  voiceActive?: boolean
  processingError?: string | null
  isModeSwitching?: boolean
  modeSwitchMessage?: string
  liveTranscript?: string
  conversationMode?: boolean
}

export function VoiceProcessingPanel({
  currentResponse,
  isThinking,
  isTyping,
  mode,
  userTranscript,
  voiceActive,
  processingError,
  isModeSwitching = false,
  modeSwitchMessage = "",
  liveTranscript = "",
  conversationMode = false
}: VoiceProcessingPanelProps) {
  const [audioWaves, setAudioWaves] = useState([1, 1, 1, 1, 1, 1, 1, 1, 1])
  const [animatedProcessingSteps, setAnimatedProcessingSteps] = useState<string[]>([])

  useEffect(() => {
    const waveInterval = setInterval(() => {
      if (isTyping || isThinking || voiceActive) {
        setAudioWaves((prev) => prev.map(() => Math.random() * 5 + 0.5))
      } else {
        setAudioWaves([1, 1, 1, 1, 1, 1, 1, 1, 1])
      }
    }, 100)
    return () => clearInterval(waveInterval)
  }, [isTyping, isThinking, voiceActive])

  useEffect(() => {
    if (isThinking) {
      const modeSpecificSteps = {
        basic: [
          "Initializing Sonar Basic search...",
          "Analyzing your query...",
          "Searching web sources...",
          "Compiling quick answers...",
          "Finalizing response..."
        ],
        reasoning: [
          "Initializing Sonar Reasoning Pro...",
          "Deep analytical processing...",
          "Applying logical reasoning...",
          "Cross-referencing knowledge base...",
          "Synthesizing thoughtful response..."
        ],
        "deep-research": [
          "Initializing Sonar Deep Research...",
          "Comprehensive source analysis...",
          "Multi-dimensional investigation...",
          "Academic source validation...",
          "Compiling research findings..."
        ]
      };
      
      const steps = modeSpecificSteps[mode];
      let stepIndex = 0
      setAnimatedProcessingSteps([steps[0]])
      const stepInterval = setInterval(() => {
        stepIndex++
        if (stepIndex < steps.length) {
          setAnimatedProcessingSteps((prev) => [...prev, steps[stepIndex]])
        } else {
          clearInterval(stepInterval)
        }
      }, 800)
      return () => clearInterval(stepInterval)
    } else {
      setAnimatedProcessingSteps([])
    }
  }, [isThinking, mode])

  const getModeConfig = () => {
    switch (mode) {
      case "reasoning":
        return {
          icon: <Brain className="w-5 h-5 text-purple-400" />,
          label: "SONAR REASONING PRO",
          gradient: "from-purple-500/20 to-purple-600/20",
          triggerHints: "Say: 'think about', 'analyze', 'figure out', 'explain why'"
        }
      case "deep-research":
        return {
          icon: <Activity className="w-5 h-5 text-blue-400" />,
          label: "SONAR DEEP RESEARCH",
          gradient: "from-blue-500/20 to-blue-600/20",
          triggerHints: "Say: 'research into', 'investigate', 'deep dive', 'study'"
        }
      default:
        return {
          icon: <Zap className="w-5 h-5 text-green-400" />,
          label: "SONAR BASIC",
          gradient: "from-green-500/20 to-green-600/20",
          triggerHints: "Say: 'search for', 'what is', 'find', 'define'"
        }
    }
  }

  const modeConfig = getModeConfig()
  const displayStatus = processingError ? "Error" : 
                       isModeSwitching ? "Switching Mode" :
                       isThinking ? "Processing" : 
                       isTyping ? "Speaking" : 
                       voiceActive ? "Listening" : "Ready"

  return (
    <div className="h-full flex flex-col bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] border border-[#20808D]/30 rounded-3xl overflow-hidden shadow-2xl backdrop-blur-xl text-[#FBFAF4]">
      {/* Header */}
      <div className={`bg-gradient-to-r ${modeConfig.gradient} border-b border-[#20808D]/30 p-6`}>
        <div className="flex items-center justify-between mb-4">
          <div className="flex items-center gap-4">
            <PerplexityLogo className="w-8 h-8" />
            <div>
              <span className="text-[#FBFAF4] font-bold text-xl">Monday</span>
              <div className="text-[#20808D] text-sm font-medium">Learning Companion</div>
            </div>
          </div>
          <div className="flex items-center gap-3 bg-[#20808D]/20 rounded-full px-4 py-2 backdrop-blur-sm">
            {modeConfig.icon}
            <span className="text-[#20808D] text-xs font-bold">{modeConfig.label}</span>
          </div>
        </div>

        {/* Enhanced Audio Visualization */}
        <div className="flex items-center justify-center gap-1 h-16 bg-[#091717]/30 rounded-2xl p-4 backdrop-blur-sm">
          {audioWaves.map((height, index) => (
            <div
              key={index}
              className="w-1.5 bg-gradient-to-t from-[#20808D] to-[#20808D]/50 rounded-full transition-all duration-150 shadow-lg"
              style={{ height: `${height * 12}px` }}
            ></div>
          ))}
        </div>
      </div>

      {/* Content Area */}
      <div className="flex-1 p-6 overflow-hidden flex flex-col">
        {/* Mode Switch Message */}
        {isModeSwitching && modeSwitchMessage && (
          <div className="mb-4 p-4 bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-400/30 rounded-xl shadow-lg">
            <div className="flex items-center gap-3">
              <div className="w-4 h-4 border-2 border-blue-400 border-t-transparent rounded-full animate-spin"></div>
              <span className="text-blue-400 font-semibold text-sm">{modeSwitchMessage}</span>
            </div>
          </div>
        )}

        {isThinking ? (
          <div className="space-y-6 flex-1 flex flex-col justify-center items-center">
            <div className="flex items-center gap-4 mb-8">
              <div className="relative">
                <Loader2 className="w-8 h-8 text-[#20808D] animate-spin" />
                {mode === 'reasoning' && <Brain className="w-4 h-4 text-purple-400 absolute -top-1 -right-1" />}
                {mode === 'deep-research' && <Activity className="w-4 h-4 text-blue-400 absolute -top-1 -right-1" />}
                {mode === 'basic' && <Zap className="w-4 h-4 text-green-400 absolute -top-1 -right-1" />}
              </div>
              <div className="flex flex-col">
                <span className="text-[#20808D] font-bold text-xl">PROCESSING YOUR REQUEST</span>
                <span className="text-[#20808D]/70 text-sm">{modeConfig.label} MODE</span>
              </div>
            </div>
            <div className="space-y-3 w-full max-w-md">
              {animatedProcessingSteps.map((step, index) => (
                <div key={index} className="flex items-center gap-3 p-2 bg-[#20808D]/10 rounded-lg animate-fade-in">
                  <div className="w-2.5 h-2.5 bg-[#20808D] rounded-full animate-pulse shadow-lg shadow-[#20808D]/50"></div>
                  <span className="text-[#FBFAF4]/90 text-sm">{step}</span>
                </div>
              ))}
            </div>
          </div>
        ) : (
          <div className="h-full flex flex-col flex-1">
            {userTranscript && (
              <div className="mb-4 p-4 bg-[#20808D]/10 border border-[#20808D]/20 rounded-xl shadow">
                <div className="flex items-center gap-3 mb-2">
                  <MessageSquare className="w-5 h-5 text-[#20808D]" />
                  <span className="text-[#20808D] font-semibold text-sm">Your Query ({modeConfig.label}):</span>
                </div>
                <p className="text-[#FBFAF4]/80 text-sm font-mono">{userTranscript}</p>
              </div>
            )}
            
            {/* Voice Trigger Hints */}
            {voiceActive && !userTranscript && !isThinking && !isTyping && (
              <div className="mb-4 p-4 bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 border border-[#20808D]/20 rounded-xl">
                <div className="flex items-center gap-3 mb-2">
                  <Mic className="w-5 h-5 text-[#20808D] animate-pulse" />
                  <span className="text-[#20808D] font-semibold text-sm">Voice Triggers for {modeConfig.label}:</span>
                </div>
                <p className="text-[#FBFAF4]/70 text-xs">{modeConfig.triggerHints}</p>
                <p className="text-[#20808D]/60 text-xs mt-1">Or just speak naturally - Monday will auto-detect the right mode!</p>
              </div>
            )}

            {/* Live Transcript Display */}
            {voiceActive && liveTranscript && !isThinking && conversationMode && (
              <div className="mb-4 p-3 bg-gradient-to-r from-green-500/20 to-green-600/20 border border-green-500/40 rounded-xl shadow-lg shadow-green-500/20">
                <div className="flex items-center gap-2 mb-2">
                  <div className="relative">
                    <div className="w-3 h-3 bg-green-400 rounded-full animate-pulse"></div>
                    <div className="absolute -inset-1 bg-green-400 rounded-full opacity-20 animate-ping"></div>
                  </div>
                  <span className="text-green-400 font-bold text-xs">LIVE RECOGNITION:</span>
                </div>
                <p className="text-green-300 text-sm font-mono font-semibold">"{liveTranscript}"</p>
                <div className="text-xs text-green-400/70 mt-1 flex items-center gap-2">
                  <div className="w-2 h-2 bg-green-400 rounded-full animate-pulse"></div>
                  Say your command...
                </div>
              </div>
            )}

            {/* Listening for Hey Monday indicator */}
            {voiceActive && !conversationMode && !isThinking && !userTranscript && (
              <div className="mb-4 p-4 bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 border border-[#20808D]/20 rounded-xl">
                <div className="flex items-center gap-3 mb-2">
                  <Mic className="w-5 h-5 text-[#20808D] animate-pulse" />
                  <span className="text-[#20808D] font-semibold text-sm">Listening for activation...</span>
                </div>
                <p className="text-[#FBFAF4]/70 text-xs">Say "Hey Monday" + your question to start</p>
                <p className="text-[#20808D]/60 text-xs mt-1">Monday will auto-detect the right mode for your question!</p>
              </div>
            )}
            
            {processingError && (
              <div className="mb-4 p-4 bg-red-500/10 border border-red-500/30 rounded-xl shadow">
                <div className="flex items-center gap-3 mb-2">
                  <Zap className="w-5 h-5 text-red-400" />
                  <span className="text-red-400 font-semibold text-sm">Error:</span>
                </div>
                <p className="text-red-300/90 text-sm font-mono">{processingError}</p>
              </div>
            )}
            <div className="flex items-center gap-4 mb-2">
              <Volume2 className="w-6 h-6 text-[#20808D]" />
              <span className="text-[#20808D] font-bold text-lg">MONDAY'S RESPONSE</span>
              {currentResponse && (
                <span className="text-[#20808D]/60 text-xs bg-[#20808D]/10 px-2 py-1 rounded-full">
                  {modeConfig.label}
                </span>
              )}
            </div>
            <div className="flex-1 bg-gradient-to-br from-[#091717] to-[#0A1A1A] rounded-2xl p-6 overflow-y-auto shadow-inner border border-[#20808D]/20 custom-scrollbar min-h-[150px]">
              {currentResponse || (!isTyping && !isThinking && !processingError) ? (
                <div className="text-[#FBFAF4] text-sm leading-relaxed whitespace-pre-wrap font-mono tracking-wide">
                  {currentResponse}
                  {isTyping && <span className="inline-block w-0.5 h-5 bg-[#20808D] ml-1 animate-pulse shadow-lg"></span>}
                </div>
              ) : (
                <div className="text-center text-[#20808D]/70 text-sm">
                  {!processingError && (voiceActive ? 
                    "🎤 Listening... Say 'Hey Monday' + your question" : 
                    "Click the mic button to start voice interaction"
                  )}
                </div>
              )}
            </div>
          </div>
        )}
      </div>

      {/* Enhanced Status Bar */}
      <div className="bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 border-t border-[#20808D]/30 p-4">
        <div className="flex items-center justify-between text-xs">
          <div className="flex items-center gap-6">
            <div className="flex items-center gap-3">
              <div
                className={`w-3 h-3 rounded-full shadow-lg ${displayStatus === "Error" ? "bg-red-500 animate-pulse shadow-red-500/50" : (isThinking || isTyping || voiceActive) ? "bg-green-400 animate-pulse shadow-green-400/50" : "bg-[#20808D] shadow-[#20808D]/50"}`}
              ></div>
              <span className="text-[#FBFAF4]/80 font-medium">{displayStatus}</span>
            </div>
            <div className="flex items-center gap-3">
              <Mic className={`w-4 h-4 ${voiceActive ? "text-green-400 animate-pulse" : "text-[#20808D]"}`} />
              <span className={`font-medium ${voiceActive ? "text-green-400" : "text-[#FBFAF4]/80"}`}>Voice {voiceActive ? "Active" : "Inactive"}</span>
            </div>
          </div>
          <div className="text-[#20808D]/80 font-medium">Powered by Perplexity</div>
        </div>
      </div>
    </div>
  )
}
