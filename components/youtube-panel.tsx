"use client"

import { useState, useEffect } from "react"
import { Play, Pause, Volume2, Maximize, ExternalLink, ThumbsUp, Share, Search } from "lucide-react"

interface YouTubePanelProps {
  videoId?: string;
  query?: string;
  mode?: "basic" | "reasoning" | "deep-research";
  videoData?: {
    title?: string;
    channel?: string;
    description?: string;
    thumbnail?: string;
    publishedAt?: string;
    relatedVideos?: Array<{
      id: string;
      title: string;
      channel: string;
      description: string;
      thumbnail: string;
      publishedAt: string;
    }>;
  };
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

export function YouTubePanel({ videoId, query, mode, videoData: backendVideoData }: YouTubePanelProps) {
  const [isPlaying, setIsPlaying] = useState(false)
  const [currentTime, setCurrentTime] = useState("0:00")
  const [duration, setDuration] = useState("0:00")
  const [videoData, setVideoData] = useState<VideoData | null>(null)
  const [relatedVideos, setRelatedVideos] = useState<VideoData[]>([])

  // Default educational video for when no data is available
  const defaultVideoData: VideoData = {
    id: "aircAruvnKk",
    title: "But what is a Neural Network? | Deep learning, chapter 1",
    description: "Subscribe for new videos every Friday! Neural networks can seem like a black box. But what exactly are they? And how do they work?",
    channelName: "3Blue1Brown",
    thumbnail: `https://img.youtube.com/vi/aircAruvnKk/maxresdefault.jpg`,
    duration: "19:13",
    publishedAt: "Oct 5, 2017",
    viewCount: "15M views"
  }

  // Effect to handle video data from backend
  useEffect(() => {
    if (videoId && backendVideoData) {
      // Use real video data from backend YouTube API
      console.log("Using real video data from YouTube API:", videoId);
      setVideoData({
        id: videoId,
        title: backendVideoData.title || `Educational Content: ${query || 'Video'}`,
        description: backendVideoData.description || `Relevant educational video content about ${query || 'this topic'}.`,
        channelName: backendVideoData.channel || "Educational Channel",
        thumbnail: backendVideoData.thumbnail || `https://img.youtube.com/vi/${videoId}/maxresdefault.jpg`,
        duration: "Educational content",
        publishedAt: backendVideoData.publishedAt || "Recently selected",
        viewCount: "Educational content"
      });
      
      // Set related videos from backend (limit to 2)
      if (backendVideoData.relatedVideos && backendVideoData.relatedVideos.length > 0) {
        const formattedRelatedVideos = backendVideoData.relatedVideos.slice(0, 2).map(video => ({
          id: video.id,
          title: video.title,
          description: video.description,
          channelName: video.channel,
          thumbnail: video.thumbnail,
          duration: "Educational content",
          publishedAt: video.publishedAt,
          viewCount: "Educational content"
        }));
        setRelatedVideos(formattedRelatedVideos);
      } else {
        setRelatedVideos([]);
      }
      
    } else if (videoId) {
      // Basic video data when we only have an ID
      console.log("Using basic video data for ID:", videoId);
      setVideoData({
        id: videoId,
        title: query ? `Educational Content: ${query}` : "Educational Content",
        description: query ? 
          `Relevant educational video content about ${query}. This video was automatically selected based on your query.` :
          "Relevant educational video content.",
        channelName: "Educational Channel",
        thumbnail: `https://img.youtube.com/vi/${videoId}/maxresdefault.jpg`,
        duration: "Educational content",
        publishedAt: "Recently selected",
        viewCount: "Educational content"
      });
      setRelatedVideos([]);
    } else {
      // Use default video when no data available
      console.log("Using default educational video");
      setVideoData(defaultVideoData);
      setRelatedVideos([]);
    }
  }, [videoId, query, backendVideoData]);

  const currentVideoId = videoData?.id || defaultVideoData.id;
  const embedUrl = `https://www.youtube.com/embed/${currentVideoId}?autoplay=0&controls=1&modestbranding=1&rel=0&cc_load_policy=1`;

  const handleVideoSelect = (video: VideoData) => {
    setVideoData(video);
    setCurrentTime("0:00");
  };

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
              <span className="text-[#FBFAF4] font-bold text-lg">
                {query ? "Educational Videos" : "Learning Content"}
              </span>
              <div className="text-[#20808D] text-sm font-medium">
                {query ? `Results for: ${query}` : "Ready for your query"}
              </div>
            </div>
          </div>
          <div className="flex items-center gap-2">
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
        {videoData ? (
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

        {videoData && (
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
        {relatedVideos.length > 0 && (
          <div className="mt-4 pt-4 border-t border-[#20808D]/20">
            <h4 className="text-[#FBFAF4] font-semibold text-sm mb-3">Related Educational Content</h4>
            <div className="space-y-2 max-h-32 overflow-y-auto">
              {relatedVideos.map((video, index) => (
                <button
                  key={video.id}
                  onClick={() => handleVideoSelect(video)}
                  className="w-full text-left p-2 rounded-lg bg-[#20808D]/10 hover:bg-[#20808D]/20 transition-colors"
                >
                  <div className="text-[#FBFAF4] text-xs font-medium truncate">{video.title}</div>
                  <div className="text-[#20808D] text-xs">{video.channelName} • {video.publishedAt}</div>
                </button>
              ))}
            </div>
          </div>
        )}
      </div>
    </div>
  )
}
