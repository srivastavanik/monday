"use client"

import { useState } from "react"
import { Play, Pause, Volume2, Maximize, ExternalLink, ThumbsUp, Share } from "lucide-react"

export function YouTubePanel() {
  const [isPlaying, setIsPlaying] = useState(false)
  const [currentTime, setCurrentTime] = useState("3:42")
  const [duration] = useState("12:18")

  return (
    <div className="h-full flex flex-col bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] border border-[#20808D]/30 rounded-3xl overflow-hidden shadow-2xl backdrop-blur-xl">
      {/* Header */}
      <div className="bg-gradient-to-r from-red-500/20 to-[#20808D]/20 border-b border-[#20808D]/30 p-6">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4">
            <div className="w-8 h-8 bg-gradient-to-br from-red-600 to-red-500 rounded-lg flex items-center justify-center shadow-lg">
              <Play className="w-4 h-4 text-white fill-white" />
            </div>
            <div>
              <span className="text-[#FBFAF4] font-bold text-lg">YouTube</span>
              <div className="text-[#20808D] text-sm font-medium">Educational Content</div>
            </div>
          </div>
          <ExternalLink className="w-5 h-5 text-[#20808D] cursor-pointer hover:text-[#FBFAF4] transition-colors" />
        </div>
      </div>

      {/* Video Area */}
      <div className="flex-1 relative bg-black rounded-t-2xl overflow-hidden">
        <iframe
          className="w-full h-full"
          src="https://www.youtube.com/embed/pYT9F8_LFTM?si=dQw4w9WgXcQ"
          title="Binary Search Algorithm Explained"
          frameBorder="0"
          allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
          allowFullScreen
        ></iframe>

        {/* Enhanced Video overlay controls */}
        <div className="absolute bottom-0 left-0 right-0 bg-gradient-to-t from-black/90 via-black/50 to-transparent p-6">
          <div className="space-y-4">
            {/* Progress bar */}
            <div className="w-full h-2 bg-white/20 rounded-full overflow-hidden backdrop-blur-sm">
              <div className="w-1/3 h-full bg-gradient-to-r from-red-600 to-red-500 rounded-full shadow-lg"></div>
            </div>

            {/* Controls */}
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-4">
                <button
                  onClick={() => setIsPlaying(!isPlaying)}
                  className="w-10 h-10 bg-white/20 hover:bg-white/30 rounded-full flex items-center justify-center transition-all backdrop-blur-sm shadow-lg"
                >
                  {isPlaying ? (
                    <Pause className="w-5 h-5 text-white" />
                  ) : (
                    <Play className="w-5 h-5 text-white ml-0.5" />
                  )}
                </button>
                <Volume2 className="w-5 h-5 text-white" />
                <span className="text-white text-sm font-medium">
                  {currentTime} / {duration}
                </span>
              </div>
              <div className="flex items-center gap-3">
                <ThumbsUp className="w-4 h-4 text-white/80 cursor-pointer hover:text-white transition-colors" />
                <Share className="w-4 h-4 text-white/80 cursor-pointer hover:text-white transition-colors" />
                <Maximize className="w-4 h-4 text-white/80 cursor-pointer hover:text-white transition-colors" />
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Enhanced Video Info */}
      <div className="bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 border-t border-[#20808D]/30 p-6">
        <h3 className="text-[#FBFAF4] font-bold text-base mb-3">Binary Search Algorithm - Complete Visual Tutorial</h3>
        <p className="text-[#FBFAF4]/80 text-sm leading-relaxed mb-4">
          Master binary search with step-by-step visualization and real-world examples. Learn the O(log n) algorithm
          that efficiently searches sorted arrays and forms the foundation of many advanced data structures.
        </p>
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs text-[#20808D] font-medium">
            <span>CS Dojo</span>
            <span>•</span>
            <span>1.2M views</span>
            <span>•</span>
            <span>2 years ago</span>
          </div>
          <div className="flex items-center gap-2">
            <div className="w-2 h-2 bg-green-400 rounded-full animate-pulse"></div>
            <span className="text-[#20808D] text-xs font-medium">Recommended</span>
          </div>
        </div>
      </div>
    </div>
  )
}
