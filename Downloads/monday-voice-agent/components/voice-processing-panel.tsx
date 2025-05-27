"use client"

import { useState, useEffect } from "react"
import { Mic, Volume2, Brain, Zap, Activity, Sparkles } from "lucide-react"
import { PerplexityLogo } from "./perplexity-logo"

interface VoiceProcessingPanelProps {
  currentResponse: string
  isThinking: boolean
  isTyping: boolean
  mode: "basic" | "reasoning" | "deep-research"
}

export function VoiceProcessingPanel({ currentResponse, isThinking, isTyping, mode }: VoiceProcessingPanelProps) {
  const [audioWaves, setAudioWaves] = useState([1, 1, 1, 1, 1, 1, 1, 1, 1])
  const [processingSteps, setProcessingSteps] = useState<string[]>([])

  useEffect(() => {
    const waveInterval = setInterval(() => {
      if (isTyping || isThinking) {
        setAudioWaves((prev) => prev.map(() => Math.random() * 5 + 0.5))
      } else {
        setAudioWaves([1, 1, 1, 1, 1, 1, 1, 1, 1])
      }
    }, 100)
    return () => clearInterval(waveInterval)
  }, [isTyping, isThinking])

  useEffect(() => {
    if (isThinking) {
      const steps = [
        "Initializing Sonar API connection...",
        "Parsing natural language query...",
        "Accessing knowledge database...",
        "Accessing knowledge database...",
        "Processing with neural networks...",
        "Generating contextual response...",
        "Optimizing for learning outcomes...",
      ]

      let stepIndex = 0
      setProcessingSteps([steps[0]])

      const stepInterval = setInterval(() => {
        stepIndex++
        if (stepIndex < steps.length) {
          setProcessingSteps((prev) => [...prev, steps[stepIndex]])
        } else {
          clearInterval(stepInterval)
        }
      }, 500)

      return () => clearInterval(stepInterval)
    } else {
      setProcessingSteps([])
    }
  }, [isThinking])

  const getModeConfig = () => {
    switch (mode) {
      case "reasoning":
        return {
          icon: <Brain className="w-5 h-5 text-[#20808D]" />,
          label: "SONAR REASONING PRO",
          gradient: "from-purple-500/20 to-[#20808D]/20",
        }
      case "deep-research":
        return {
          icon: <Activity className="w-5 h-5 text-[#20808D]" />,
          label: "SONAR DEEP RESEARCH",
          gradient: "from-blue-500/20 to-[#20808D]/20",
        }
      default:
        return {
          icon: <Zap className="w-5 h-5 text-[#20808D]" />,
          label: "SONAR BASIC",
          gradient: "from-green-500/20 to-[#20808D]/20",
        }
    }
  }

  const modeConfig = getModeConfig()

  return (
    <div className="h-full flex flex-col bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] border border-[#20808D]/30 rounded-3xl overflow-hidden shadow-2xl backdrop-blur-xl">
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
      <div className="flex-1 p-6 overflow-hidden">
        {isThinking ? (
          <div className="space-y-6">
            <div className="flex items-center gap-4 mb-8">
              <div className="relative">
                <Brain className="w-8 h-8 text-[#20808D] animate-pulse" />
                <Sparkles className="w-4 h-4 text-[#20808D] absolute -top-1 -right-1 animate-bounce" />
              </div>
              <span className="text-[#20808D] font-bold text-xl">PROCESSING</span>
            </div>

            <div className="space-y-4">
              {processingSteps.map((step, index) => (
                <div key={index} className="flex items-center gap-4 animate-fade-in">
                  <div className="w-3 h-3 bg-[#20808D] rounded-full animate-pulse shadow-lg shadow-[#20808D]/50"></div>
                  <span className="text-[#FBFAF4]/90 text-sm font-medium">{step}</span>
                </div>
              ))}
            </div>

            {/* Enhanced thinking animation */}
            <div className="mt-12 flex justify-center">
              <div className="relative">
                <div className="w-20 h-20 border-2 border-[#20808D]/30 rounded-full"></div>
                <div className="absolute inset-0 w-20 h-20 border-2 border-[#20808D] rounded-full border-t-transparent animate-spin"></div>
                <div className="absolute inset-2 w-16 h-16 border border-[#20808D]/50 rounded-full border-b-transparent animate-spin animate-reverse"></div>
                <div className="absolute inset-4 w-12 h-12 border border-[#20808D]/30 rounded-full border-l-transparent animate-spin"></div>
              </div>
            </div>
          </div>
        ) : (
          <div className="h-full flex flex-col">
            <div className="flex items-center gap-4 mb-6">
              <Volume2 className="w-6 h-6 text-[#20808D]" />
              <span className="text-[#20808D] font-bold text-lg">RESPONSE</span>
            </div>

            <div className="flex-1 bg-gradient-to-br from-[#091717] to-[#0A1A1A] rounded-2xl p-6 overflow-y-auto shadow-inner border border-[#20808D]/20 custom-scrollbar">
              <div className="text-[#20808D] text-sm leading-relaxed whitespace-pre-wrap font-mono tracking-wide">
                {currentResponse || ""}
                {isTyping && <span className="inline-block w-0.5 h-5 bg-[#20808D] ml-1 animate-pulse shadow-lg"></span>}
              </div>
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
                className={`w-3 h-3 rounded-full shadow-lg ${
                  isThinking || isTyping
                    ? "bg-green-400 animate-pulse shadow-green-400/50"
                    : "bg-[#20808D] shadow-[#20808D]/50"
                }`}
              ></div>
              <span className="text-[#FBFAF4]/80 font-medium">
                {isThinking ? "Processing" : isTyping ? "Speaking" : "Ready"}
              </span>
            </div>
            <div className="flex items-center gap-3">
              <Mic className="w-4 h-4 text-[#20808D]" />
              <span className="text-[#FBFAF4]/80 font-medium">Voice Active</span>
            </div>
          </div>
          <div className="text-[#20808D]/80 font-medium">Powered by Perplexity</div>
        </div>
      </div>
    </div>
  )
}
