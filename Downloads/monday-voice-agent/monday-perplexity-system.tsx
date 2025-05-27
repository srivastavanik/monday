"use client"

import { useState, useEffect } from "react"
import { BinaryTree3D } from "./components/binary-tree-3d"
import { VoiceProcessingPanel } from "./components/voice-processing-panel"
import { YouTubePanel } from "./components/youtube-panel"
import { PerplexityLogo } from "./components/perplexity-logo"
import { Brain, Search, Activity, Zap, Sparkles } from "lucide-react"

const responses = {
  basic: `Binary search is a fundamental algorithm that efficiently locates elements in sorted arrays through a divide-and-conquer approach. The algorithm works by repeatedly comparing the target value with the middle element of the current search space.

Here's how it operates:
1. Start with the entire sorted array
2. Calculate the middle index
3. Compare target with middle element
4. If equal, return the index
5. If target is smaller, search left half
6. If target is larger, search right half
7. Repeat until found or search space is empty

Time Complexity: O(log n) - exponentially faster than linear search
Space Complexity: O(1) iterative, O(log n) recursive

The beauty of binary search lies in its logarithmic efficiency. In an array of 1 million elements, it requires at most 20 comparisons to find any element, compared to 500,000 average comparisons for linear search.

This algorithm forms the foundation for many advanced data structures including binary search trees, B-trees used in databases, and various optimization algorithms in machine learning.`,

  reasoning: `Let me analyze why binary search is particularly effective for the binary search tree structure you're observing in the visualization.

The fundamental insight is the ordering property: every node maintains the invariant that left children are smaller and right children are larger. This creates a natural search path that binary search can exploit.

Consider searching for value 60:
1. Start at root (50) - target is larger, so the answer must be in the right subtree
2. Move to node (70) - target is smaller, so it must be in the left subtree  
3. Reach node (60) - exact match found!

Each comparison eliminates approximately half the remaining nodes, achieving O(log n) performance. The tree structure makes this elimination visual and intuitive - we can literally see which branches become irrelevant.

This spatial representation helps understand why binary search works: it's not just about array indices, but about systematically narrowing down possibilities using the ordering property. The algorithm's efficiency comes from this guaranteed reduction in search space at each step.

The visualization demonstrates how abstract algorithmic concepts become concrete when represented spatially, making the learning process more intuitive and memorable.`,

  research: `Binary search represents a cornerstone achievement in algorithmic design, embodying the divide-and-conquer paradigm that has influenced countless computational solutions across computer science.

Historical Development:
The concept emerged from John Mauchly's 1946 work, but the first fully correct implementation wasn't published until 1962 by Donald Knuth. This 16-year gap illustrates the subtle complexity in handling edge cases - a reminder that even "simple" algorithms require careful implementation.

Mathematical Foundation:
The algorithm's O(log n) complexity stems from the recurrence relation T(n) = T(n/2) + O(1). Each iteration halves the problem size, leading to logarithmic depth. This mathematical property makes binary search optimal for comparison-based searching in sorted data.

Variants and Extensions:
- Lower bound search: Finding the first occurrence of a value
- Upper bound search: Finding the insertion point for maintaining sorted order
- Exponential search: Handling unbounded arrays
- Interpolation search: Exploiting uniform distribution for O(log log n) performance
- Binary search on answer: Solving optimization problems by searching the solution space

Real-world Applications:
Database systems use binary search in B-tree indices for rapid record retrieval. Git's bisect command employs binary search to locate bug-introducing commits. Operating systems use it for memory allocation and process scheduling. Machine learning algorithms apply binary search for hyperparameter optimization and feature selection.

The binary search tree visualization demonstrates the algorithm's spatial nature - how logical comparisons translate to physical navigation through data structures. This spatial understanding is crucial for grasping more complex tree-based algorithms like AVL trees, red-black trees, and B-trees used in modern databases.

Contemporary Relevance:
In our era of big data, binary search's logarithmic scaling becomes even more valuable. While linear algorithms become impractical with massive datasets, binary search maintains efficiency, making it indispensable for real-time systems and large-scale applications.`,
}

export default function MondayPerplexitySystem() {
  const [currentMode, setCurrentMode] = useState<"basic" | "reasoning" | "deep-research">("basic")
  const [currentResponse, setCurrentResponse] = useState("")
  const [fullResponse, setFullResponse] = useState(responses.basic || "")
  const [isThinking, setIsThinking] = useState(false)
  const [isTyping, setIsTyping] = useState(false)
  const [charIndex, setCharIndex] = useState(0)

  // Typewriter effect
  useEffect(() => {
    if (isTyping && fullResponse && charIndex < fullResponse.length) {
      const timeout = setTimeout(
        () => {
          setCurrentResponse(fullResponse.slice(0, charIndex + 1))
          setCharIndex(charIndex + 1)
        },
        25 + Math.random() * 15,
      )
      return () => clearTimeout(timeout)
    } else if (fullResponse && charIndex >= fullResponse.length) {
      setIsTyping(false)
    }
  }, [charIndex, fullResponse, isTyping])

  const startNewResponse = (mode: "basic" | "reasoning" | "deep-research") => {
    const responseText = responses[mode] || ""
    setCurrentMode(mode)
    setFullResponse(responseText)
    setIsThinking(true)
    setCurrentResponse("")
    setCharIndex(0)

    setTimeout(() => {
      setIsThinking(false)
      setIsTyping(true)
    }, 3500)
  }

  useEffect(() => {
    startNewResponse("basic")
  }, [])

  const getModeConfig = (mode: string) => {
    switch (mode) {
      case "reasoning":
        return {
          icon: <Brain className="w-4 h-4" />,
          label: "Reasoning Pro",
          gradient: "from-purple-500/30 to-[#20808D]/30",
        }
      case "deep-research":
        return {
          icon: <Activity className="w-4 h-4" />,
          label: "Deep Research",
          gradient: "from-blue-500/30 to-[#20808D]/30",
        }
      default:
        return {
          icon: <Search className="w-4 h-4" />,
          label: "Basic",
          gradient: "from-green-500/30 to-[#20808D]/30",
        }
    }
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] p-8">
      <div className="w-full max-w-[2000px] mx-auto">
        {/* Enhanced Header */}
        <div className="text-center mb-12">
          <div className="flex items-center justify-center gap-6 mb-6">
            <PerplexityLogo className="w-16 h-16" />
            <div className="text-center">
              <h1 className="text-5xl font-bold text-[#FBFAF4] tracking-tight mb-2">Monday</h1>
              <p className="text-[#20808D] text-lg font-semibold">Learning, Amalgamated</p>
            </div>
            <div className="flex items-center gap-3 bg-gradient-to-r from-[#20808D]/20 to-[#20808D]/10 rounded-2xl px-6 py-3 backdrop-blur-xl border border-[#20808D]/30">
              <Zap className="w-5 h-5 text-[#20808D]" />
              <span className="text-[#20808D] text-sm font-bold">VR AI Companion</span>
              <Sparkles className="w-4 h-4 text-[#20808D] animate-pulse" />
            </div>
          </div>
        </div>

        {/* Enhanced Mode Selector */}
        <div className="flex justify-center gap-4 mb-12">
          {(["basic", "reasoning", "deep-research"] as const).map((mode) => {
            const config = getModeConfig(mode)
            return (
              <button
                key={mode}
                onClick={() => startNewResponse(mode)}
                className={`flex items-center gap-3 px-6 py-3 rounded-2xl transition-all text-sm font-bold backdrop-blur-xl border ${
                  currentMode === mode
                    ? "bg-gradient-to-r from-[#20808D] to-[#20808D]/80 text-[#FBFAF4] border-[#20808D]/50 shadow-lg shadow-[#20808D]/30"
                    : "bg-gradient-to-r from-[#20808D]/20 to-[#20808D]/10 text-[#20808D] hover:from-[#20808D]/30 hover:to-[#20808D]/20 border-[#20808D]/30"
                }`}
              >
                {config.icon}
                Sonar {config.label}
              </button>
            )
          })}
        </div>

        {/* Main Content Grid with Concave Effect */}
        <div className="grid grid-cols-1 xl:grid-cols-12 gap-8 h-[800px]">
          {/* Left Panel - 3D Visualization (Angled inward) */}
          <div className="xl:col-span-4 transform xl:-rotate-2 xl:scale-95 origin-right">
            <div className="h-full shadow-2xl shadow-[#20808D]/20">
              <BinaryTree3D />
            </div>
          </div>

          {/* Center Panel - Voice Processing (Normal) */}
          <div className="xl:col-span-4 z-10">
            <div className="h-full shadow-2xl shadow-[#20808D]/30">
              <VoiceProcessingPanel
                currentResponse={currentResponse}
                isThinking={isThinking}
                isTyping={isTyping}
                mode={currentMode}
              />
            </div>
          </div>

          {/* Right Panel - YouTube (Angled inward) */}
          <div className="xl:col-span-4 transform xl:rotate-2 xl:scale-95 origin-left">
            <div className="h-full shadow-2xl shadow-[#20808D]/20">
              <YouTubePanel />
            </div>
          </div>
        </div>

        {/* Enhanced Footer */}
        <div className="text-center mt-12">
          <div className="inline-flex items-center gap-4 bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 rounded-2xl px-8 py-4 backdrop-blur-xl border border-[#20808D]/20">
            <PerplexityLogo className="w-6 h-6" />
            <span className="text-[#20808D]/80 text-sm font-medium">
              Powered by Perplexity Sonar API • WebXR Ready • Quest 2/3 Compatible
            </span>
          </div>
        </div>
      </div>
    </div>
  )
}
