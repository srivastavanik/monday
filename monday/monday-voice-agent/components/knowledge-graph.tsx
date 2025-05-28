"use client"

import type { Source } from "../types/monday"
import { ExternalLink, Star } from "lucide-react"

interface KnowledgeGraphProps {
  sources: Source[]
  query: string
}

export function KnowledgeGraph({ sources, query }: KnowledgeGraphProps) {
  return (
    <div className="space-y-6">
      <div className="text-center mb-8">
        <h3 className="text-blue-300 font-semibold mb-2">RESEARCH CONSTELLATION</h3>
        <p className="text-slate-400 text-sm">Query: "{query}"</p>
      </div>

      {/* Central query node */}
      <div className="relative flex justify-center mb-8">
        <div className="w-32 h-32 bg-gradient-to-br from-blue-500 to-blue-600 rounded-full flex items-center justify-center border-4 border-blue-400 shadow-lg shadow-blue-500/30">
          <span className="text-white font-bold text-center text-sm px-2">
            {query.length > 30 ? query.substring(0, 30) + "..." : query}
          </span>
        </div>
      </div>

      {/* Source nodes */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        {sources.map((source, index) => (
          <div key={source.id} className="relative">
            {/* Connection line to center */}
            <div className="absolute top-4 left-1/2 w-px h-8 bg-blue-400/30 transform -translate-x-1/2 -translate-y-full"></div>

            <div className="bg-slate-800/60 border border-blue-500/20 rounded-lg p-4 hover:border-blue-400/40 transition-all duration-300">
              <div className="flex items-start justify-between mb-2">
                <h4 className="text-white font-semibold text-sm leading-tight flex-1">{source.title}</h4>
                <ExternalLink className="w-4 h-4 text-blue-400 ml-2 flex-shrink-0" />
              </div>

              <p className="text-slate-300 text-xs mb-3 leading-relaxed">{source.excerpt}</p>

              {/* Relevance indicator */}
              <div className="flex items-center justify-between">
                <div className="flex items-center gap-1">
                  <Star className="w-3 h-3 text-yellow-400" />
                  <span className="text-xs text-slate-400">{Math.round(source.relevance * 100)}% relevant</span>
                </div>

                <div className="w-16 h-1 bg-slate-700 rounded-full overflow-hidden">
                  <div
                    className="h-full bg-gradient-to-r from-yellow-400 to-green-400 transition-all duration-1000"
                    style={{ width: `${source.relevance * 100}%` }}
                  ></div>
                </div>
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
