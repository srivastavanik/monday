"use client"

import type { LearningNode } from "../types/monday"
import { Brain, Zap } from "lucide-react"

interface LearningConstellationProps {
  nodes: LearningNode[]
  activeNode?: string
}

export function LearningConstellation({ nodes, activeNode }: LearningConstellationProps) {
  return (
    <div className="relative w-full h-64 bg-slate-900/30 rounded-lg border border-blue-500/20 overflow-hidden">
      <div className="absolute inset-0 bg-gradient-to-br from-blue-950/20 to-transparent"></div>

      <div className="absolute top-4 left-4 flex items-center gap-2">
        <Brain className="w-4 h-4 text-blue-400" />
        <span className="text-blue-300 text-sm font-semibold">LEARNING CONSTELLATION</span>
      </div>

      <div className="relative w-full h-full p-8">
        {nodes.map((node) => {
          const isActive = node.id === activeNode
          const opacity = isActive ? 1 : 0.3 + node.understanding * 0.7
          const size = isActive ? 16 : 8 + node.understanding * 8

          return (
            <div
              key={node.id}
              className="absolute transform -translate-x-1/2 -translate-y-1/2 transition-all duration-500"
              style={{
                left: `${node.position.x}%`,
                top: `${node.position.y}%`,
                opacity,
              }}
            >
              {/* Node */}
              <div
                className={`rounded-full border-2 flex items-center justify-center transition-all duration-500 ${
                  isActive
                    ? "bg-blue-500 border-blue-400 shadow-lg shadow-blue-500/50"
                    : "bg-blue-600/60 border-blue-500/60"
                }`}
                style={{ width: size, height: size }}
              >
                {isActive && <Zap className="w-3 h-3 text-white" />}
              </div>

              {/* Label */}
              <div className="absolute top-full left-1/2 transform -translate-x-1/2 mt-1">
                <span className="text-xs text-blue-300 whitespace-nowrap bg-slate-900/80 px-2 py-1 rounded">
                  {node.topic}
                </span>
              </div>

              {/* Understanding indicator */}
              <div className="absolute -bottom-1 left-1/2 transform -translate-x-1/2">
                <div className="w-8 h-1 bg-slate-700 rounded-full overflow-hidden">
                  <div
                    className="h-full bg-gradient-to-r from-blue-400 to-green-400 transition-all duration-1000"
                    style={{ width: `${node.understanding * 100}%` }}
                  ></div>
                </div>
              </div>
            </div>
          )
        })}

        {/* Connection lines */}
        <svg className="absolute inset-0 w-full h-full pointer-events-none">
          {nodes.map((node) =>
            node.connections.map((connectionId) => {
              const connectedNode = nodes.find((n) => n.id === connectionId)
              if (!connectedNode) return null

              return (
                <line
                  key={`${node.id}-${connectionId}`}
                  x1={`${node.position.x}%`}
                  y1={`${node.position.y}%`}
                  x2={`${connectedNode.position.x}%`}
                  y2={`${connectedNode.position.y}%`}
                  stroke="rgb(59 130 246 / 0.3)"
                  strokeWidth="1"
                  className="transition-all duration-500"
                />
              )
            }),
          )}
        </svg>
      </div>
    </div>
  )
}
