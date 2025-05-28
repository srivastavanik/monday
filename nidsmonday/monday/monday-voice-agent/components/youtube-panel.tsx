"use client"

import { useState, useEffect } from "react"
import { Play, Pause, Volume2, Maximize, ExternalLink, ThumbsUp, Share, Search, RefreshCw } from "lucide-react"

interface YouTubePanelProps {
  videoId?: string;
  query?: string;
  mode?: "basic" | "reasoning" | "deep-research";
}

interface VideoData {
  id: string;
  title: string;
  description: string;
  channelName: string;
  thumbnail: string;
  duration: string;
  publishedAt: string;
  viewCount: string;
}

export function YouTubePanel({ videoId, query, mode }: YouTubePanelProps) {
  const [isPlaying, setIsPlaying] = useState(false)
  const [currentTime, setCurrentTime] = useState("0:00")
  const [duration, setDuration] = useState("0:00")
  const [videoData, setVideoData] = useState<VideoData | null>(null)
  const [isSearching, setIsSearching] = useState(false)
  const [searchError, setSearchError] = useState<string | null>(null)
  const [relatedVideos, setRelatedVideos] = useState<VideoData[]>([])

  // Default educational video for demo
  const defaultVideoData: VideoData = {
    id: "pYT9F8_LFTM",
    title: "Binary Search Algorithm - Complete Visual Tutorial",
    description: "Master binary search with step-by-step visualization and real-world examples. Learn the O(log n) algorithm that efficiently searches sorted arrays and forms the foundation of many advanced data structures.",
    channelName: "CS Dojo",
    thumbnail: `https://img.youtube.com/vi/pYT9F8_LFTM/maxresdefault.jpg`,
    duration: "12:18",
    publishedAt: "2 years ago",
    viewCount: "1.2M views"
  }

  // Search for educational videos based on query
  const searchEducationalVideos = async (searchQuery: string) => {
    if (!searchQuery) return;
    
    setIsSearching(true);
    setSearchError(null);
    
    try {
      // For demo purposes, we'll simulate video search with educational content
      // In production, you would use YouTube Data API v3
      const educationalQueries = [
        `${searchQuery} tutorial`,
        `${searchQuery} explained`,
        `${searchQuery} course`,
        `${searchQuery} lecture`,
        `learn ${searchQuery}`,
        `${searchQuery} fundamentals`
      ];
      
      // Simulate API delay
      await new Promise(resolve => setTimeout(resolve, 1000));
      
      // Generate mock educational video data based on query
      const mockVideos: VideoData[] = educationalQueries.map((q, index) => ({
        id: `mock_${Date.now()}_${index}`,
        title: `${searchQuery.charAt(0).toUpperCase() + searchQuery.slice(1)} - ${
          ['Complete Tutorial', 'Explained Simply', 'Full Course', 'Comprehensive Guide', 'Step by Step', 'Fundamentals'][index]
        }`,
        description: `Learn ${searchQuery} with this comprehensive educational video. Perfect for ${
          mode === 'reasoning' ? 'deep analytical thinking' : 
          mode === 'deep-research' ? 'thorough research and understanding' : 
          'foundational learning'
        }.`,
        channelName: ['Khan Academy', 'MIT OpenCourseWare', 'Coursera', 'edX', 'FreeCodeCamp', 'Crash Course'][index],
        thumbnail: `https://img.youtube.com/vi/dQw4w9WgXcQ/maxresdefault.jpg`,
        duration: ['15:42', '23:18', '45:30', '12:05', '38:22', '19:47'][index],
        publishedAt: ['1 year ago', '6 months ago', '2 years ago', '3 months ago', '1 month ago', '8 months ago'][index],
        viewCount: ['2.1M', '856K', '3.4M', '1.7M', '945K', '1.3M'][index] + ' views'
      }));
      
      setRelatedVideos(mockVideos);
      
      // Set the first video as the main video
      if (mockVideos.length > 0) {
        setVideoData(mockVideos[0]);
      }
      
    } catch (error) {
      console.error('Error searching for videos:', error);
      setSearchError('Failed to search for educational videos');
      setVideoData(defaultVideoData);
    } finally {
      setIsSearching(false);
    }
  };

  // Effect to search for videos when query changes
  useEffect(() => {
    if (query && query.trim()) {
      searchEducationalVideos(query.trim());
    } else if (videoId) {
      // If a specific video ID is provided, use it
      setVideoData({
        id: videoId,
        title: "Educational Content",
        description: "Relevant educational video content.",
        channelName: "Educational Channel",
        thumbnail: `https://img.youtube.com/vi/${videoId}/maxresdefault.jpg`,
        duration: "Unknown",
        publishedAt: "Recently",
        viewCount: "Educational content"
      });
    } else {
      // Use default video
      setVideoData(defaultVideoData);
    }
  }, [query, videoId]);

  const currentVideoId = videoData?.id || defaultVideoData.id;
  const embedUrl = `https://www.youtube.com/embed/${currentVideoId}?autoplay=0&controls=1&modestbranding=1&rel=0&cc_load_policy=1`;

  const handleVideoSelect = (video: VideoData) => {
    setVideoData(video);
    setCurrentTime("0:00");
  };

  const handleRefreshSearch = () => {
    if (query) {
      searchEducationalVideos(query);
    }
  };

  return (
    <div className="h-full flex flex-col bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] border border-[#20808D]/30 rounded-3xl overflow-hidden shadow-2xl backdrop-blur-xl">
      {/* Header */}
      <div className="bg-gradient-to-r from-red-500/20 to-[#20808D]/20 border-b border-[#20808D]/30 p-6">
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4">
            <div className="w-8 h-8 bg-gradient-to-br from-red-600 to-red-500 rounded-lg flex items-center justify-center shadow-lg">
              {isSearching ? (
                <RefreshCw className="w-4 h-4 text-white animate-spin" />
              ) : (
                <Play className="w-4 h-4 text-white fill-white" />
              )}
            </div>
            <div>
              <span className="text-[#FBFAF4] font-bold text-lg">
                {query ? "Educational Videos" : "Learning Content"}
              </span>
              <div className="text-[#20808D] text-sm font-medium">
                {query ? `Results for: ${query}` : "Ready for your query"}
              </div>
            </div>
          </div>
          <div className="flex items-center gap-2">
            {query && (
              <button
                onClick={handleRefreshSearch}
                disabled={isSearching}
                className="p-2 rounded-lg bg-[#20808D]/20 hover:bg-[#20808D]/30 transition-colors"
                title="Refresh search"
              >
                <RefreshCw className={`w-4 h-4 text-[#20808D] ${isSearching ? 'animate-spin' : ''}`} />
              </button>
            )}
            {videoData && (
              <a 
                href={`https://www.youtube.com/watch?v=${currentVideoId}`} 
                target="_blank" 
                rel="noopener noreferrer" 
                title="Open in YouTube"
                className="p-2 rounded-lg bg-[#20808D]/20 hover:bg-[#20808D]/30 transition-colors"
              >
                <ExternalLink className="w-4 h-4 text-[#20808D]" />
              </a>
            )}
          </div>
        </div>
      </div>

      {/* Video Area */}
      <div className="flex-1 relative bg-black rounded-t-2xl overflow-hidden">
        {isSearching ? (
          <div className="w-full h-full flex items-center justify-center bg-gradient-to-br from-gray-900 to-black">
            <div className="text-center">
              <RefreshCw className="w-8 h-8 text-[#20808D] animate-spin mx-auto mb-4" />
              <p className="text-[#20808D] text-sm">Searching for educational videos...</p>
            </div>
          </div>
        ) : searchError ? (
          <div className="w-full h-full flex items-center justify-center bg-gradient-to-br from-red-900/20 to-black">
            <div className="text-center">
              <Search className="w-8 h-8 text-red-400 mx-auto mb-4" />
              <p className="text-red-400 text-sm">{searchError}</p>
              <button
                onClick={handleRefreshSearch}
                className="mt-2 px-4 py-2 bg-red-500/20 hover:bg-red-500/30 rounded-lg text-red-400 text-xs transition-colors"
              >
                Try Again
              </button>
            </div>
          </div>
        ) : videoData ? (
          <iframe
            className="w-full h-full"
            src={embedUrl}
            title={videoData.title}
            frameBorder="0"
            allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
            allowFullScreen
          ></iframe>
        ) : (
          <div className="w-full h-full flex items-center justify-center bg-gradient-to-br from-gray-900 to-black">
            <div className="text-center">
              <Search className="w-8 h-8 text-gray-500 mx-auto mb-4" />
              <p className="text-gray-500 text-sm">Ask Monday a question to find educational videos</p>
            </div>
          </div>
        )}

        {videoData && !isSearching && (
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
                    {currentTime} / {videoData.duration}
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
        )}
      </div>

      {/* Video Info */}
      <div className="bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 border-t border-[#20808D]/30 p-6">
        {videoData && (
          <>
            <h3 className="text-[#FBFAF4] font-bold text-base mb-3">{videoData.title}</h3>
            <p className="text-[#FBFAF4]/80 text-sm leading-relaxed mb-4 line-clamp-3">
              {videoData.description}
            </p>
            <div className="flex items-center justify-between">
              <div className="flex items-center gap-4 text-xs text-[#20808D] font-medium">
                <span>{videoData.channelName}</span>
                <span>•</span>
                <span>{videoData.viewCount}</span>
                <span>•</span>
                <span>{videoData.publishedAt}</span>
              </div>
              <div className="flex items-center gap-2">
                <div className="w-2 h-2 bg-green-400 rounded-full animate-pulse"></div>
                <span className="text-[#20808D] text-xs font-medium">Educational</span>
              </div>
            </div>
          </>
        )}
        
        {/* Related Videos */}
        {relatedVideos.length > 1 && (
          <div className="mt-4 pt-4 border-t border-[#20808D]/20">
            <h4 className="text-[#FBFAF4] font-semibold text-sm mb-3">Related Educational Content</h4>
            <div className="space-y-2 max-h-32 overflow-y-auto">
              {relatedVideos.slice(1, 4).map((video, index) => (
                <button
                  key={video.id}
                  onClick={() => handleVideoSelect(video)}
                  className="w-full text-left p-2 rounded-lg bg-[#20808D]/10 hover:bg-[#20808D]/20 transition-colors"
                >
                  <div className="text-[#FBFAF4] text-xs font-medium truncate">{video.title}</div>
                  <div className="text-[#20808D] text-xs">{video.channelName} • {video.duration}</div>
                </button>
              ))}
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
