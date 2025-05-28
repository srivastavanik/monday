"use client"

import { useState, useEffect } from "react"
import { Play, Pause, Volume2, Maximize, ExternalLink, ThumbsUp, Share, Search, RefreshCw, Youtube } from "lucide-react"

interface YouTubePanelProps {
  videoId?: string
  query?: string
  mode?: "basic" | "reasoning" | "deep-research"
  responseContent?: string  // Add the actual model response content
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

// YouTube API configuration
const YOUTUBE_API_KEY = "AIzaSyDn_zCV8AGjkQFufH7RDGkSiXD75-2Q39M";
const YOUTUBE_API_BASE_URL = "https://www.googleapis.com/youtube/v3";

// Placeholder component for YouTube panel
function YouTubePlaceholder() {
  return (
    <div className="w-full h-full flex items-center justify-center bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] rounded-2xl border border-[#20808D]/30">
      <div className="text-center">
        <div className="mb-4">
          <div className="w-20 h-20 mx-auto bg-[#20808D]/20 rounded-full flex items-center justify-center">
            <Youtube className="w-10 h-10 text-[#20808D]" />
          </div>
        </div>
        <h3 className="text-lg font-semibold text-[#20808D] mb-2">Query Monday to find relevant videos</h3>
        <p className="text-sm text-[#20808D]/60 max-w-xs mx-auto">
          I'll search for educational videos that complement your learning journey
        </p>
      </div>
    </div>
  )
}

export function YouTubePanel({ videoId, query, mode, responseContent }: YouTubePanelProps) {
  // All hooks must be declared at the top, before any conditional returns
  const [videoData, setVideoData] = useState<VideoData | null>(null)
  const [isSearching, setIsSearching] = useState(false)
  const [searchError, setSearchError] = useState<string | null>(null)
  const [relatedVideos, setRelatedVideos] = useState<VideoData[]>([])
  const [selectedVideoId, setSelectedVideoId] = useState<string | null>(null)

  // Function to format ISO 8601 duration to human readable
  const formatDuration = (isoDuration: string): string => {
    try {
      const match = isoDuration.match(/PT(?:(\d+)H)?(?:(\d+)M)?(?:(\d+)S)?/);
      if (!match) return "Unknown";
      
      const hours = match[1] ? parseInt(match[1]) : 0;
      const minutes = match[2] ? parseInt(match[2]) : 0;
      const seconds = match[3] ? parseInt(match[3]) : 0;
      
      if (hours > 0) {
        return `${hours}:${minutes.toString().padStart(2, '0')}:${seconds.toString().padStart(2, '0')}`;
      }
      return `${minutes}:${seconds.toString().padStart(2, '0')}`;
    } catch (error) {
      return "Unknown";
    }
  };

  // Function to format view count
  const formatViewCount = (count: string): string => {
    const num = parseInt(count);
    if (num >= 1000000) {
      return `${(num / 1000000).toFixed(1)}M views`;
    } else if (num >= 1000) {
      return `${(num / 1000).toFixed(0)}K views`;
    }
    return `${num} views`;
  };

  // Function to calculate relative time
  const getRelativeTime = (dateString: string): string => {
    try {
      const date = new Date(dateString);
      const now = new Date();
      const diffMs = now.getTime() - date.getTime();
      const diffDays = Math.floor(diffMs / (1000 * 60 * 60 * 24));
      
      if (diffDays === 0) return "Today";
      if (diffDays === 1) return "Yesterday";
      if (diffDays < 7) return `${diffDays} days ago`;
      if (diffDays < 30) return `${Math.floor(diffDays / 7)} weeks ago`;
      if (diffDays < 365) return `${Math.floor(diffDays / 30)} months ago`;
      return `${Math.floor(diffDays / 365)} years ago`;
    } catch (error) {
      return "Recently";
    }
  };

  // Function to extract keywords from model response content
  const extractKeywords = (content: string, originalQuery: string): string[] => {
    if (!content) return [originalQuery];
    
    // Remove common words and extract meaningful terms
    const stopWords = new Set([
      'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for', 'of', 'with', 'by',
      'is', 'are', 'was', 'were', 'be', 'been', 'being', 'have', 'has', 'had', 'do', 'does', 'did',
      'will', 'would', 'could', 'should', 'may', 'might', 'can', 'this', 'that', 'these', 'those',
      'it', 'its', 'they', 'them', 'their', 'we', 'us', 'our', 'you', 'your', 'i', 'me', 'my'
    ]);
    
    // Extract words from content, prioritizing nouns and important terms
    const words = content.toLowerCase()
      .replace(/[^\w\s]/g, ' ')
      .split(/\s+/)
      .filter(word => word.length > 2 && !stopWords.has(word));
    
    // Get word frequency
    const wordFreq: { [key: string]: number } = {};
    words.forEach(word => {
      wordFreq[word] = (wordFreq[word] || 0) + 1;
    });
    
    // Sort by frequency and take top keywords
    const keywords = Object.entries(wordFreq)
      .sort(([,a], [,b]) => b - a)
      .slice(0, 5)
      .map(([word]) => word);
    
    // Always include the original query as the primary keyword
    const finalKeywords = [originalQuery, ...keywords].filter((keyword, index, arr) => 
      arr.indexOf(keyword) === index
    );
    
    return finalKeywords.slice(0, 6); // Limit to 6 keywords max
  };

  // Function to build enhanced search query from keywords
  const buildSearchQuery = (keywords: string[], mode: string): string => {
    const primaryKeyword = keywords[0] || '';
    const additionalKeywords = keywords.slice(1, 3); // Use 2 additional keywords
    
    let baseQuery = primaryKeyword;
    if (additionalKeywords.length > 0) {
      baseQuery += ' ' + additionalKeywords.join(' ');
    }
    
    // Add mode-specific terms
    switch (mode) {
      case 'reasoning':
        return `${baseQuery} explained tutorial analysis how why`;
      case 'deep-research':
        return `${baseQuery} documentary comprehensive research study`;
      default:
        return `${baseQuery} tutorial guide explanation`;
    }
  };

  // Search for videos using YouTube Data API
  const searchYouTubeVideos = async (searchQuery: string, modelResponse?: string) => {
    if (!searchQuery) return;
    
    setIsSearching(true);
    setSearchError(null);
    
    try {
      console.log('Searching for videos with query:', searchQuery);
      console.log('Model response available:', !!modelResponse);
      
      // Extract keywords from model response if available
      const keywords = extractKeywords(modelResponse || '', searchQuery);
      const enhancedQuery = buildSearchQuery(keywords, mode || 'basic');
      
      console.log('Extracted keywords:', keywords);
      console.log('Enhanced search query:', enhancedQuery);
      
      // Search for videos with enhanced query
      let searchData;
      let searchResponse = await fetch(
        `${YOUTUBE_API_BASE_URL}/search?` + new URLSearchParams({
          part: 'snippet',
          q: enhancedQuery,
          type: 'video',
          maxResults: '10', // Get more to filter better
          order: 'relevance',
          videoDuration: 'medium', // Stick to medium length videos
          videoDefinition: 'any',
          videoLicense: 'any',
          key: YOUTUBE_API_KEY
        })
      );
      
      if (!searchResponse.ok) {
        const errorData = await searchResponse.text();
        console.error('YouTube search API error:', errorData);
        
        // Try a simpler fallback search if the enhanced query failed
        if (enhancedQuery !== searchQuery) {
          console.log('Trying fallback search with original query...');
          searchResponse = await fetch(
            `${YOUTUBE_API_BASE_URL}/search?` + new URLSearchParams({
              part: 'snippet',
              q: searchQuery + ' tutorial',
              type: 'video',
              maxResults: '5',
              order: 'relevance',
              key: YOUTUBE_API_KEY
            })
          );
          
          if (!searchResponse.ok) {
            throw new Error(`Failed to search videos: ${searchResponse.status} - ${errorData}`);
          }
        } else {
          throw new Error(`Failed to search videos: ${searchResponse.status} - ${errorData}`);
        }
      }
      
      searchData = await searchResponse.json();
      
      if (!searchData.items || searchData.items.length === 0) {
        console.log('No videos found for query:', enhancedQuery);
        throw new Error('No videos found for this topic');
      }
      
      console.log(`Found ${searchData.items.length} videos`);
      
      // Get video IDs for details
      const videoIds = searchData.items.map((item: any) => item.id.videoId).join(',');
      
      // Get detailed video information to validate videos exist
      const detailsResponse = await fetch(
        `${YOUTUBE_API_BASE_URL}/videos?` + new URLSearchParams({
          part: 'snippet,contentDetails,statistics',
          id: videoIds,
          key: YOUTUBE_API_KEY
        })
      );
      
      if (!detailsResponse.ok) {
        throw new Error('Failed to get video details');
      }
      
      const detailsData = await detailsResponse.json();
      console.log(`Got details for ${detailsData.items.length} videos`);
      
      // Filter for valid videos with simpler validation
      const validVideos: VideoData[] = detailsData.items
        .filter((item: any) => {
          const snippet = item.snippet;
          const statistics = item.statistics;
          const contentDetails = item.contentDetails;
          
          // Basic validation - just check if video exists and has basic data
          const hasValidTitle = snippet && snippet.title && !snippet.title.toLowerCase().includes('deleted');
          const hasValidChannel = snippet && snippet.channelTitle && snippet.channelTitle !== '[Deleted Channel]';
          const hasViews = statistics && statistics.viewCount && parseInt(statistics.viewCount) > 0;
          const hasDuration = contentDetails && contentDetails.duration && contentDetails.duration !== 'PT0S';
          
          console.log(`Video ${item.id} validation:`, {
            hasValidTitle,
            hasValidChannel,
            hasViews,
            hasDuration,
            title: snippet?.title
          });
          
          return hasValidTitle && hasValidChannel && hasViews && hasDuration;
        })
        .map((item: any) => ({
          id: item.id,
          title: item.snippet.title,
          description: item.snippet.description || 'Educational video content',
          channelName: item.snippet.channelTitle,
          thumbnail: item.snippet.thumbnails.medium?.url || item.snippet.thumbnails.default?.url,
          duration: formatDuration(item.contentDetails.duration),
          publishedAt: getRelativeTime(item.snippet.publishedAt),
          viewCount: formatViewCount(item.statistics.viewCount)
        }))
        .slice(0, 3); // Limit to exactly 3 videos
      
      if (validVideos.length === 0) {
        throw new Error('No valid videos found. Please try a different search term.');
      }
      
      console.log(`Filtered to ${validVideos.length} valid videos:`, validVideos.map(v => v.title));
      
      setRelatedVideos(validVideos);
      
      // Set the first video as selected
      if (validVideos.length > 0) {
        console.log('Setting first video as selected:', validVideos[0].title);
        setSelectedVideoId(validVideos[0].id);
        setVideoData(validVideos[0]);
      }
      
    } catch (error) {
      console.error('Error searching YouTube videos:', error);
      setSearchError(error instanceof Error ? error.message : 'Failed to search for videos. Please try again.');
    } finally {
      setIsSearching(false);
    }
  };

  // Effect to search for videos when query changes
  useEffect(() => {
    if (query && query.trim()) {
      // Use both the query and response content for better keyword extraction
      searchYouTubeVideos(query.trim(), responseContent);
    } else if (videoId) {
      // If a specific video ID is provided, use it
      setSelectedVideoId(videoId);
      
      // Create a temporary video data object
      const tempVideoData = {
        id: videoId,
        title: 'Loading video information...',
        description: '',
        channelName: '',
        thumbnail: '',
        duration: '',
        publishedAt: '',
        viewCount: ''
      };
      
      setVideoData(tempVideoData);
      
      // Try to get video details from YouTube API
      fetch(`${YOUTUBE_API_BASE_URL}/videos?part=snippet,contentDetails,statistics&id=${videoId}&key=${YOUTUBE_API_KEY}`)
        .then(response => response.json())
        .then(data => {
          if (data.items && data.items.length > 0) {
            const videoDetails = data.items[0];
            const formattedData = {
              id: videoDetails.id,
              title: videoDetails.snippet.title,
              description: videoDetails.snippet.description,
              channelName: videoDetails.snippet.channelTitle,
              thumbnail: videoDetails.snippet.thumbnails.high?.url || videoDetails.snippet.thumbnails.default.url,
              duration: formatDuration(videoDetails.contentDetails.duration),
              publishedAt: getRelativeTime(videoDetails.snippet.publishedAt),
              viewCount: formatViewCount(videoDetails.statistics.viewCount)
            };
            setVideoData(formattedData);
          }
        })
        .catch(error => {
          console.error('Error fetching video details:', error);
        });
    }
  }, [query, videoId, mode, responseContent]); // Add responseContent as dependency

  const getModeLabel = () => {
    switch(mode) {
      case "reasoning": return "Reasoning Mode"
      case "deep-research": return "Deep Research"
      default: return "Basic Search"
    }
  };

  const handleVideoSelect = (video: VideoData) => {
    setVideoData(video);
    setSelectedVideoId(video.id);
  };

  const handleRefreshSearch = () => {
    if (query) {
      searchYouTubeVideos(query, responseContent);
    }
  };

  // If no video ID or query yet, show placeholder
  if (!videoId && !query) {
    return <YouTubePlaceholder />
  }

  const currentVideoId = selectedVideoId || videoId;

  return (
    <div className="w-full h-full bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] rounded-2xl border border-[#20808D]/30 overflow-hidden">
      <div className="h-full flex flex-col">
        {/* Header */}
        <div className="bg-gradient-to-b from-[#091717]/90 to-transparent p-4 border-b border-[#20808D]/20">
          <div className="flex items-center justify-between mb-2">
            <div className="flex items-center gap-2">
              <Youtube className="w-5 h-5 text-[#20808D]" />
              <h3 className="text-[#20808D] font-bold text-sm">Educational Videos</h3>
            </div>
            <div className="flex items-center gap-2">
              <span className="text-xs text-[#20808D]/60">{getModeLabel()}</span>
              {query && (
                <button
                  onClick={handleRefreshSearch}
                  disabled={isSearching}
                  className="p-1.5 bg-[#20808D]/20 hover:bg-[#20808D]/30 rounded-lg transition-colors disabled:opacity-50"
                  title="Refresh search"
                >
                  <RefreshCw className={`w-3 h-3 text-[#FBFAF4] ${isSearching ? 'animate-spin' : ''}`} />
                </button>
              )}
            </div>
          </div>
          {query && (
            <p className="text-[#FBFAF4]/60 text-xs truncate">
              <Search className="w-3 h-3 inline mr-1" />
              {query}
            </p>
          )}
        </div>

        {/* Video Player Area */}
        <div className="flex-1 relative">
          {searchError ? (
            <div className="w-full h-full flex items-center justify-center p-6">
              <div className="text-center">
                <div className="mb-4 p-3 bg-red-500/10 rounded-full inline-block">
                  <Youtube className="w-8 h-8 text-red-400" />
                </div>
                <p className="text-red-400 text-sm mb-3">{searchError}</p>
                {query && (
                  <button
                    onClick={handleRefreshSearch}
                    className="px-4 py-2 bg-[#20808D]/20 hover:bg-[#20808D]/30 rounded-lg transition-colors text-[#FBFAF4] text-sm"
                  >
                    Try Again
                  </button>
                )}
              </div>
            </div>
          ) : isSearching ? (
            <div className="w-full h-full flex items-center justify-center">
              <div className="text-center">
                <div className="mb-4">
                  <div className="w-12 h-12 border-3 border-[#20808D] border-t-transparent rounded-full animate-spin mx-auto"></div>
                </div>
                <p className="text-[#20808D] text-sm">Searching for educational videos...</p>
              </div>
            </div>
          ) : currentVideoId ? (
            <div className="w-full h-full bg-black rounded-lg overflow-hidden">
              <iframe
                src={`https://www.youtube.com/embed/${currentVideoId}?autoplay=0&controls=1&modestbranding=1&rel=0&fs=1`}
                title={videoData?.title || 'Educational Video'}
                className="w-full h-full"
                frameBorder="0"
                allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share"
                allowFullScreen
                onError={(e) => {
                  console.error('YouTube iframe error for video:', currentVideoId);
                  setSearchError('This video cannot be played. Trying next video...');
                  
                  // Try next video if current one fails
                  const currentIndex = relatedVideos.findIndex(v => v.id === currentVideoId);
                  if (currentIndex >= 0 && currentIndex < relatedVideos.length - 1) {
                    setTimeout(() => {
                      const nextVideo = relatedVideos[currentIndex + 1];
                      handleVideoSelect(nextVideo);
                    }, 2000);
                  }
                }}
              />
              
              {/* Video Info Overlay */}
              {videoData && (
                <div className="absolute bottom-0 left-0 right-0 bg-gradient-to-t from-black/90 via-black/50 to-transparent p-4">
                  <div className="flex items-center justify-between">
                    <div className="flex-1 min-w-0">
                      <h3 className="text-white font-semibold text-sm truncate">{videoData.title}</h3>
                      <p className="text-white/80 text-xs">{videoData.channelName} • {videoData.viewCount}</p>
                    </div>
                    <a
                      href={`https://www.youtube.com/watch?v=${currentVideoId}`}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="ml-3 px-3 py-1 bg-red-600 hover:bg-red-700 text-white text-xs rounded transition-colors"
                    >
                      View on YouTube
                    </a>
                  </div>
                </div>
              )}
            </div>
          ) : (
            <YouTubePlaceholder />
          )}
        </div>
        
        {/* Related Videos - Show exactly 2 additional videos (excluding current) */}
        {relatedVideos.length > 1 && (
          <div className="border-t border-[#20808D]/20 p-3 bg-[#091717]/80">
            <h4 className="text-[#20808D] text-sm font-semibold mb-3">Related Videos</h4>
            <div className="space-y-2">
              {relatedVideos
                .filter(video => video.id !== currentVideoId) // Exclude currently playing video
                .slice(0, 2) // Show only 2 additional videos
                .map((video) => (
                  <div
                    key={video.id}
                    className="flex gap-3 p-2 rounded-lg bg-[#20808D]/10 hover:bg-[#20808D]/20 border border-[#20808D]/30 hover:border-[#20808D]/50 transition-all cursor-pointer"
                    onClick={() => handleVideoSelect(video)}
                  >
                    <img
                      src={video.thumbnail}
                      alt={video.title}
                      className="w-20 h-12 object-cover rounded flex-shrink-0"
                      onError={(e) => {
                        (e.target as HTMLImageElement).src = `https://img.youtube.com/vi/${video.id}/mqdefault.jpg`;
                      }}
                    />
                    <div className="flex-1 min-w-0">
                      <p className="text-[#FBFAF4] text-sm font-medium line-clamp-2 mb-1">{video.title}</p>
                      <p className="text-[#20808D]/70 text-xs">{video.channelName}</p>
                      <div className="flex items-center justify-between mt-1">
                        <span className="text-[#20808D]/60 text-xs">{video.duration} • {video.viewCount}</span>
                        <a
                          href={`https://www.youtube.com/watch?v=${video.id}`}
                          target="_blank"
                          rel="noopener noreferrer"
                          onClick={(e) => e.stopPropagation()}
                          className="text-[#20808D] hover:text-[#FBFAF4] text-xs underline"
                        >
                          View on YouTube
                        </a>
                      </div>
                    </div>
                  </div>
                ))}
            </div>
            {relatedVideos.length === 1 && (
              <p className="text-[#20808D]/60 text-xs text-center py-4">
                Only one video found for this topic
              </p>
            )}
          </div>
        )}
      </div>
    </div>
  )
}
