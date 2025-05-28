"use client"

import { useState, useEffect } from "react"
import { Mic, Volume2, Zap, Brain } from "lucide-react"

const voiceResponses = [
  "Good morning! I'm here to help you tackle your Monday goals.",
  "I've analyzed your schedule and found 3 optimization opportunities.",
  "Your productivity metrics show a 23% improvement this week.",
  "I've prepared your daily briefing with 7 priority items.",
  "Weather looks great today - perfect for that outdoor meeting.",
  "I've detected 4 new messages that require your attention.",
  "Your calendar has been updated with the latest changes.",
  "I'm monitoring your energy levels and suggest a 10-minute break.",
  "I've found 2 relevant articles based on your recent searches.",
  "Your next meeting starts in 15 minutes. Shall I prepare the materials?",
]

export default function MondayVoiceAgent() {
  const [currentResponse, setCurrentResponse] = useState("")
  const [fullResponse, setFullResponse] = useState(voiceResponses[0])
  const [isThinking, setIsThinking] = useState(false)
  const [isTyping, setIsTyping] = useState(false)
  const [audioWaves, setAudioWaves] = useState([1, 1, 1, 1, 1])
  const [charIndex, setCharIndex] = useState(0)

  // Typewriter effect
  useEffect(() => {
    if (isTyping && charIndex < fullResponse.length) {
      const timeout = setTimeout(
        () => {
          setCurrentResponse(fullResponse.slice(0, charIndex + 1))
          setCharIndex(charIndex + 1)
        },
        50 + Math.random() * 50,
      ) // Variable typing speed for more natural feel

      return () => clearTimeout(timeout)
    } else if (charIndex >= fullResponse.length) {
      setIsTyping(false)
    }
  }, [charIndex, fullResponse, isTyping])

  // Main response cycle
  useEffect(() => {
    const startTyping = () => {
      setIsThinking(true)
      setCurrentResponse("")
      setCharIndex(0)

      // Thinking delay
      setTimeout(() => {
        setIsThinking(false)
        setIsTyping(true)
      }, 1500)
    }

    // Start first response
    startTyping()

    const interval = setInterval(() => {
      const randomIndex = Math.floor(Math.random() * voiceResponses.length)
      setFullResponse(voiceResponses[randomIndex])
      startTyping()
    }, 8000)

    return () => clearInterval(interval)
  }, [])

  // Audio waves animation
  useEffect(() => {
    const waveInterval = setInterval(() => {
      if (isTyping) {
        setAudioWaves((prev) => prev.map(() => Math.random() * 3 + 0.5))
      } else {
        setAudioWaves([1, 1, 1, 1, 1])
      }
    }, 150)

    return () => clearInterval(waveInterval)
  }, [isTyping])

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-950 via-blue-950 to-slate-900 flex items-center justify-center p-6">
      <div className="w-full max-w-4xl">
        {/* Header */}
        <div className="text-center mb-12">
          <div className="flex items-center justify-center gap-3 mb-4">
            <div className="relative">
              <div className="w-12 h-12 bg-blue-500 rounded-full flex items-center justify-center">
                <Zap className="w-6 h-6 text-white" />
              </div>
              <div className="absolute -top-1 -right-1 w-4 h-4 bg-green-400 rounded-full animate-pulse"></div>
            </div>
            <h1 className="text-4xl font-bold text-white tracking-tight">Monday</h1>
          </div>
          <p className="text-blue-300 text-lg font-light">AI Voice Assistant</p>
        </div>

        {/* Main Voice Panel */}
        <div className="relative">
          {/* Background glow */}
          <div className="absolute inset-0 bg-blue-500/10 rounded-3xl blur-xl"></div>

          <div className="relative bg-slate-900/80 backdrop-blur-sm border border-blue-500/20 rounded-3xl p-8 shadow-2xl">
            {/* Audio Visualization */}
            <div className="flex items-center justify-center gap-1 mb-8">
              {audioWaves.map((height, index) => (
                <div
                  key={index}
                  className={`w-1 rounded-full transition-all duration-150 ${
                    isTyping ? "bg-blue-400" : "bg-slate-600"
                  }`}
                  style={{ height: `${height * 20}px` }}
                ></div>
              ))}
            </div>

            {/* Voice Response */}
            <div className="text-center min-h-[120px] flex flex-col justify-center">
              <div className="flex items-center justify-center gap-3 mb-6">
                {isThinking ? (
                  <Brain className="w-5 h-5 text-blue-400 animate-pulse" />
                ) : (
                  <Volume2 className="w-5 h-5 text-blue-400" />
                )}
                <span className="text-blue-300 font-semibold text-lg tracking-wide">MONDAY:</span>
              </div>

              <div className="relative">
                {isThinking ? (
                  <div className="flex items-center justify-center gap-2">
                    <div className="flex gap-1">
                      <div className="w-2 h-2 bg-blue-400 rounded-full animate-bounce"></div>
                      <div
                        className="w-2 h-2 bg-blue-400 rounded-full animate-bounce"
                        style={{ animationDelay: "0.1s" }}
                      ></div>
                      <div
                        className="w-2 h-2 bg-blue-400 rounded-full animate-bounce"
                        style={{ animationDelay: "0.2s" }}
                      ></div>
                    </div>
                    <span className="text-slate-400 text-sm ml-2">thinking...</span>
                  </div>
                ) : (
                  <div className="text-white text-xl leading-relaxed font-light max-w-3xl mx-auto">
                    {currentResponse}
                    {isTyping && <span className="inline-block w-0.5 h-6 bg-blue-400 ml-1 animate-pulse"></span>}
                  </div>
                )}
              </div>
            </div>

            {/* Status Indicators */}
            <div className="flex items-center justify-center gap-6 mt-8">
              <div className="flex items-center gap-2">
                <div
                  className={`w-2 h-2 rounded-full ${isThinking || isTyping ? "bg-green-400 animate-pulse" : "bg-slate-600"}`}
                ></div>
                <span className="text-slate-400 text-sm">
                  {isThinking ? "Thinking" : isTyping ? "Speaking" : "Idle"}
                </span>
              </div>
              <div className="flex items-center gap-2">
                <Mic className="w-4 h-4 text-blue-400" />
                <span className="text-slate-400 text-sm">Listening</span>
              </div>
              <div className="flex items-center gap-2">
                <div className={`w-2 h-2 rounded-full ${isTyping ? "bg-blue-400" : "bg-slate-600"}`}></div>
                <span className="text-slate-400 text-sm">Voice Active</span>
              </div>
            </div>
          </div>
        </div>

        {/* Bottom Info */}
        <div className="text-center mt-8">
          <p className="text-slate-500 text-sm">
            Voice agent powered by advanced AI • Always learning, always improving
          </p>
        </div>
      </div>
    </div>
  )
}
