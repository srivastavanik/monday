"use client"

import { useState, useEffect } from "react"
import { Mic, Volume2, Zap, Brain, Search, Focus, History } from "lucide-react"
import { sampleResponses } from "./data/sample-responses"
import { ReasoningVisualization } from "./components/reasoning-visualization"
import { KnowledgeGraph } from "./components/knowledge-graph"
import { LearningConstellation } from "./components/learning-constellation"
import type { APIResponse, LearningNode } from "./types/monday"

const learningNodes: LearningNode[] = [
  { id: "1", topic: "Machine Learning", understanding: 0.8, connections: ["2", "3"], position: { x: 30, y: 40 } },
  { id: "2", topic: "Neural Networks", understanding: 0.9, connections: ["1", "4"], position: { x: 60, y: 30 } },
  { id: "3", topic: "Algorithms", understanding: 0.7, connections: ["1", "5"], position: { x: 20, y: 70 } },
  { id: "4", topic: "Deep Learning", understanding: 0.6, connections: ["2"], position: { x: 80, y: 50 } },
  { id: "5", topic: "Data Structures", understanding: 0.85, connections: ["3"], position: { x: 40, y: 80 } },
]

export default function MondayLearningSystem() {
  const [currentScreen, setCurrentScreen] = useState<"main" | "reasoning" | "research">("main")
  const [currentResponse, setCurrentResponse] = useState<APIResponse>(sampleResponses[0])
  const [displayedText, setDisplayedText] = useState("")
  const [isThinking, setIsThinking] = useState(false)
  const [isTyping, setIsTyping] = useState(false)
  const [charIndex, setCharIndex] = useState(0)
  const [currentReasoningStep, setCurrentReasoningStep] = useState(0)
  const [audioWaves, setAudioWaves] = useState([1, 1, 1, 1, 1])

  // Typewriter effect for long responses
  useEffect(() => {
    if (isTyping && charIndex < currentResponse.content.length) {
      const timeout = setTimeout(
        () => {
          setDisplayedText(currentResponse.content.slice(0, charIndex + 1))
          setCharIndex(charIndex + 1)
        },
        20 + Math.random() * 30, // Faster for long content
      )
      return () => clearTimeout(timeout)
    } else if (charIndex >= currentResponse.content.length) {
      setIsTyping(false)
    }
  }, [charIndex, currentResponse.content, isTyping])

  // Reasoning step progression
  useEffect(() => {
    if (currentScreen === "reasoning" && currentResponse.reasoning_steps) {
      const interval = setInterval(() => {
        setCurrentReasoningStep((prev) => {
          if (prev < currentResponse.reasoning_steps!.length - 1) {
            return prev + 1
          }
          return prev
        })
      }, 2000)
      return () => clearInterval(interval)
    }
  }, [currentScreen, currentResponse.reasoning_steps])

  // Audio waves animation
  useEffect(() => {
    const waveInterval = setInterval(() => {
      if (isTyping || isThinking) {
        setAudioWaves((prev) => prev.map(() => Math.random() * 3 + 0.5))
      } else {
        setAudioWaves([1, 1, 1, 1, 1])
      }
    }, 150)
    return () => clearInterval(waveInterval)
  }, [isTyping, isThinking])

  const startNewResponse = (responseIndex: number) => {
    const response = sampleResponses[responseIndex]
    setCurrentResponse(response)
    setIsThinking(true)
    setDisplayedText("")
    setCharIndex(0)
    setCurrentReasoningStep(0)

    // Set screen based on response mode
    if (response.mode === "reasoning") {
      setCurrentScreen("reasoning")
    } else if (response.mode === "deep-research") {
      setCurrentScreen("research")
    } else {
      setCurrentScreen("main")
    }

    // Thinking delay
    setTimeout(() => {
      setIsThinking(false)
      setIsTyping(true)
    }, 2000)
  }

  const renderMainScreen = () => (
    <div className="space-y-8">
      {/* Voice Response Panel */}
      <div className="relative">
        <div className="absolute inset-0 bg-blue-500/10 rounded-3xl blur-xl"></div>
        <div className="relative bg-slate-900/80 backdrop-blur-sm border border-blue-500/20 rounded-3xl p-8 shadow-2xl">
          {/* Audio Visualization */}
          <div className="flex items-center justify-center gap-1 mb-8">
            {audioWaves.map((height, index) => (
              <div
                key={index}
                className={`w-1 rounded-full transition-all duration-150 ${
                  isTyping || isThinking ? "bg-blue-400" : "bg-slate-600"
                }`}
                style={{ height: `${height * 20}px` }}
              ></div>
            ))}
          </div>

          {/* Response Content */}
          <div className="text-center min-h-[200px] flex flex-col justify-center">
            <div className="flex items-center justify-center gap-3 mb-6">
              {isThinking ? (
                <Brain className="w-5 h-5 text-blue-400 animate-pulse" />
              ) : (
                <Volume2 className="w-5 h-5 text-blue-400" />
              )}
              <span className="text-blue-300 font-semibold text-lg tracking-wide">MONDAY:</span>
            </div>

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
                <span className="text-slate-400 text-sm ml-2">processing your query...</span>
              </div>
            ) : (
              <div className="text-white text-lg leading-relaxed font-light max-w-4xl mx-auto text-left">
                <div className="bg-slate-800/40 rounded-lg p-6 max-h-96 overflow-y-auto">
                  <div className="whitespace-pre-wrap">
                    {displayedText}
                    {isTyping && <span className="inline-block w-0.5 h-6 bg-blue-400 ml-1 animate-pulse"></span>}
                  </div>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>

      {/* Learning Constellation */}
      <LearningConstellation nodes={learningNodes} activeNode={currentResponse.mode === "basic" ? "1" : undefined} />
    </div>
  )

  const renderReasoningScreen = () => (
    <div className="space-y-8">
      {/* Reasoning Visualization */}
      <div className="relative">
        <div className="absolute inset-0 bg-blue-500/10 rounded-3xl blur-xl"></div>
        <div className="relative bg-slate-900/80 backdrop-blur-sm border border-blue-500/20 rounded-3xl p-8 shadow-2xl">
          {currentResponse.reasoning_steps && (
            <ReasoningVisualization steps={currentResponse.reasoning_steps} currentStep={currentReasoningStep} />
          )}
        </div>
      </div>

      {/* Overall Confidence */}
      <div className="text-center">
        <div className="inline-flex items-center gap-3 bg-slate-800/60 rounded-full px-6 py-3 border border-blue-500/20">
          <Brain className="w-5 h-5 text-blue-400" />
          <span className="text-white font-semibold">Overall Confidence:</span>
          <div className="w-24 h-2 bg-slate-700 rounded-full overflow-hidden">
            <div
              className="h-full bg-gradient-to-r from-blue-500 to-green-400 transition-all duration-1000"
              style={{ width: `${(currentResponse.confidence || 0) * 100}%` }}
            ></div>
          </div>
          <span className="text-blue-300 font-bold">{Math.round((currentResponse.confidence || 0) * 100)}%</span>
        </div>
      </div>
    </div>
  )

  const renderResearchScreen = () => (
    <div className="space-y-8">
      {/* Research Content */}
      <div className="relative">
        <div className="absolute inset-0 bg-blue-500/10 rounded-3xl blur-xl"></div>
        <div className="relative bg-slate-900/80 backdrop-blur-sm border border-blue-500/20 rounded-3xl p-8 shadow-2xl">
          <div className="text-white text-lg leading-relaxed font-light">
            <div className="bg-slate-800/40 rounded-lg p-6 max-h-64 overflow-y-auto mb-6">
              <div className="whitespace-pre-wrap">
                {displayedText}
                {isTyping && <span className="inline-block w-0.5 h-6 bg-blue-400 ml-1 animate-pulse"></span>}
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Knowledge Graph */}
      <div className="relative">
        <div className="absolute inset-0 bg-blue-500/5 rounded-3xl blur-xl"></div>
        <div className="relative bg-slate-900/60 backdrop-blur-sm border border-blue-500/20 rounded-3xl p-8 shadow-2xl">
          {currentResponse.sources && (
            <KnowledgeGraph sources={currentResponse.sources} query={currentResponse.query} />
          )}
        </div>
      </div>
    </div>
  )

  return (
    <div className="min-h-screen bg-gradient-to-br from-slate-950 via-blue-950 to-slate-900 p-6">
      <div className="w-full max-w-6xl mx-auto">
        {/* Header */}
        <div className="text-center mb-8">
          <div className="flex items-center justify-center gap-3 mb-4">
            <div className="relative">
              <div className="w-12 h-12 bg-blue-500 rounded-full flex items-center justify-center">
                <Zap className="w-6 h-6 text-white" />
              </div>
              <div className="absolute -top-1 -right-1 w-4 h-4 bg-green-400 rounded-full animate-pulse"></div>
            </div>
            <h1 className="text-4xl font-bold text-white tracking-tight">Monday</h1>
            <span className="text-blue-300 text-sm bg-blue-500/20 px-3 py-1 rounded-full">Learning, Amalgamated</span>
          </div>
          <p className="text-blue-300 text-lg font-light">VR AI Learning Companion</p>
        </div>

        {/* Mode Selector */}
        <div className="flex justify-center gap-4 mb-8">
          <button
            onClick={() => startNewResponse(0)}
            className={`flex items-center gap-2 px-4 py-2 rounded-lg transition-all ${
              currentScreen === "main" ? "bg-blue-500 text-white" : "bg-slate-800 text-slate-300 hover:bg-slate-700"
            }`}
          >
            <Search className="w-4 h-4" />
            Basic Query
          </button>
          <button
            onClick={() => startNewResponse(1)}
            className={`flex items-center gap-2 px-4 py-2 rounded-lg transition-all ${
              currentScreen === "reasoning"
                ? "bg-blue-500 text-white"
                : "bg-slate-800 text-slate-300 hover:bg-slate-700"
            }`}
          >
            <Brain className="w-4 h-4" />
            Reasoning Mode
          </button>
          <button
            onClick={() => startNewResponse(2)}
            className={`flex items-center gap-2 px-4 py-2 rounded-lg transition-all ${
              currentScreen === "research" ? "bg-blue-500 text-white" : "bg-slate-800 text-slate-300 hover:bg-slate-700"
            }`}
          >
            <History className="w-4 h-4" />
            Deep Research
          </button>
        </div>

        {/* Content */}
        {currentScreen === "main" && renderMainScreen()}
        {currentScreen === "reasoning" && renderReasoningScreen()}
        {currentScreen === "research" && renderResearchScreen()}

        {/* Status Bar */}
        <div className="flex items-center justify-center gap-6 mt-8">
          <div className="flex items-center gap-2">
            <div
              className={`w-2 h-2 rounded-full ${isThinking || isTyping ? "bg-green-400 animate-pulse" : "bg-slate-600"}`}
            ></div>
            <span className="text-slate-400 text-sm">
              {isThinking ? "Processing" : isTyping ? "Responding" : "Ready"}
            </span>
          </div>
          <div className="flex items-center gap-2">
            <Mic className="w-4 h-4 text-blue-400" />
            <span className="text-slate-400 text-sm">Voice Active</span>
          </div>
          <div className="flex items-center gap-2">
            <Focus className="w-4 h-4 text-blue-400" />
            <span className="text-slate-400 text-sm">VR Mode</span>
          </div>
        </div>
      </div>
    </div>
  )
}
