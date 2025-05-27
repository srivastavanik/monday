"use client"

import { useState, useEffect, useCallback } from "react"
import { BinaryTree3D } from "./components/binary-tree-3d"
import { VoiceProcessingPanel } from "./components/voice-processing-panel"
import { YouTubePanel } from "./components/youtube-panel"
import { PerplexityLogo } from "./components/perplexity-logo"
import { Brain, Search, Activity, Zap, Sparkles, Mic, MicOff } from "lucide-react"
import useMondayWebSocket from "./hooks/useMondayWebSocket"
import SpeechRecognition, { useSpeechRecognition } from "react-speech-recognition"
import { useTextToSpeech } from "./hooks/useTextToSpeech"

// Define types for API responses and WebSocket messages
interface BackendApiResponse {
  answer: string;
  youtubeVideoId?: string;
  visualizationData?: any;
  modelUsed: "basic" | "reasoning" | "deep-research";
  panels?: SpatialPanel[];
  citations?: any[];
  reasoning?: any[];
  metadata?: any;
}

interface SpatialPanel {
  id: string;
  type: string;
  position: [number, number, number];
  rotation: [number, number, number];
  title: string;
  content: string;
  fullContent?: string;
  isActive: boolean;
  opacity: number;
  createdAt: number;
  model: string;
  citations?: any[];
  reasoning?: any[];
  sources?: any[];
  isThinking?: boolean;
  isResearching?: boolean;
}

interface WebSocketMessage {
  type: "voice_response" | "reasoning_progress" | "research_progress" | "voice_error" | "spatial_response";
  message?: string;
  data?: {
    panels?: SpatialPanel[];
    mode?: string;
    model?: string;
    query?: string;
    citations?: any[];
    reasoning?: any[];
    metadata?: any;
  };
  error?: string;
  update?: any;
  query?: string;
  model?: string;
  timestamp?: number;
}

const MONDAY_WEBSOCKET_URL = "http://localhost:3001";

export default function MondayPerplexitySystem() {
  const [currentMode, setCurrentMode] = useState<"basic" | "reasoning" | "deep-research">("basic")
  const [currentResponse, setCurrentResponse] = useState("")
  const [fullResponseText, setFullResponseText] = useState("")
  const [isThinking, setIsThinking] = useState(false)
  const [isReceiving, setIsReceiving] = useState(false)
  const [charIndex, setCharIndex] = useState(0)
  const [youtubeVideoId, setYoutubeVideoId] = useState<string | undefined>(undefined)
  const [visualizationData, setVisualizationData] = useState<any>(null)
  const [userTranscript, setUserTranscript] = useState("")
  const [apiError, setApiError] = useState<string | null>(null)
  const [isMounted, setIsMounted] = useState(false)
  const [spatialPanels, setSpatialPanels] = useState<SpatialPanel[]>([])
  const [isInConversation, setIsInConversation] = useState(false)
  const [progressUpdates, setProgressUpdates] = useState<string[]>([])
  const [isProcessingVoice, setIsProcessingVoice] = useState(false)
  const [isListeningForTrigger, setIsListeningForTrigger] = useState(true)
  const [lastTriggerTime, setLastTriggerTime] = useState(0)
  const [conversationMode, setConversationMode] = useState(false)
  const [modeSwitchMessage, setModeSwitchMessage] = useState<string>("")
  const [isModeSwitching, setIsModeSwitching] = useState(false)

  // --- Text-to-Speech Hook with ElevenLabs API ---
  const { speak, stop: stopSpeaking, isSpeaking, initializeAudioContext } = useTextToSpeech({
    voiceId: "21m00Tcm4TlvDq8ikWAM", // Rachel voice
    apiKey: "sk_9454bca5ee475d45dc6e50bcb33bc3fb76f138e2191ff47d"
  });

  // --- WebSocket Hook --- 
  const { 
    isConnected: wsIsConnected,
    lastMessage: wsLastMessage,
    sendMessage: wsSendMessage,
    error: wsError 
  } = useMondayWebSocket({
    socketUrl: MONDAY_WEBSOCKET_URL,
    onOpen: () => {
      console.log("WebSocket connected to Monday backend")
      setApiError(null)
    },
    onClose: () => {
      console.log("WebSocket disconnected from Monday backend")
      setApiError("Connection to Monday backend lost. Attempting to reconnect...")
    },
    onError: (event) => {
      console.error("WebSocket error:", event)
      setApiError("Failed to connect to Monday backend. Please check if the backend is running.")
    },
  });

  // --- Speech Recognition Hook ---
  const {
    transcript,
    listening,
    resetTranscript,
    browserSupportsSpeechRecognition,
    isMicrophoneAvailable
  } = useSpeechRecognition();

  // Fix hydration issue by tracking mount state
  useEffect(() => {
    setIsMounted(true)
    // Initialize audio context and start automatic listening on mount
    if (typeof window !== 'undefined') {
      initializeAudioContext();
      // Start listening immediately after component mounts
      if (browserSupportsSpeechRecognition && isMicrophoneAvailable) {
        startAutomaticListening();
      }
    }
  }, [])

  // Monitor transcript for "Hey Monday" trigger or continuous conversation
  useEffect(() => {
    if (transcript && listening) {
      const transcriptLower = transcript.toLowerCase().trim();
      
      if (!conversationMode) {
        // Check for "Hey Monday" activation
        if (transcriptLower.includes("hey monday") || transcriptLower.includes("hey monty")) {
          const currentTime = Date.now();
          
          // Prevent duplicate triggers within 1 second (reduced from 3 seconds)
          if (currentTime - lastTriggerTime > 1000) {
            setLastTriggerTime(currentTime);
            handleMondayActivation(transcript);
          }
        }
      } else {
        // In conversation mode - send everything to Perplexity after a pause
        // Use a simple word count threshold to detect when user finishes speaking
        const wordCount = transcriptLower.split(/\s+/).filter(word => word.length > 0).length;
        
        if (wordCount >= 3) { // Wait for at least 3 words before processing
          // Set a timeout to process the command after user stops speaking (reduced from 2s to 1s)
          const timeoutId = setTimeout(() => {
            if (transcript === transcriptLower) { // Check if transcript hasn't changed
              handleContinuousConversation(transcript);
            }
          }, 1000); // 1 second pause before processing (reduced from 2000ms)
          
          return () => clearTimeout(timeoutId);
        }
      }
    }
  }, [transcript, listening, conversationMode, lastTriggerTime]);

  const startAutomaticListening = useCallback(() => {
    if (!browserSupportsSpeechRecognition || !isMicrophoneAvailable) {
      setApiError("Voice recognition not available. Please check microphone permissions.");
      return;
    }

    if (!listening) {
      SpeechRecognition.startListening({ 
        continuous: true,
        language: 'en-US'
      });
      setIsListeningForTrigger(true);
      console.log("Started automatic listening for voice input");
    }
  }, [browserSupportsSpeechRecognition, isMicrophoneAvailable, listening]);

  const handleMondayActivation = useCallback((fullTranscript: string) => {
    console.log("Hey Monday activation detected:", fullTranscript);
    
    // Stop any ongoing speech immediately
    if (isSpeaking) {
      stopSpeaking();
    }

    // Extract command after "Hey Monday"
    const triggerIndex = fullTranscript.toLowerCase().indexOf("hey monday");
    const commandPart = fullTranscript.substring(triggerIndex + 10).trim();
    
    if (commandPart && commandPart.length > 0) {
      // Hey Monday with immediate command - process instantly
      setUserTranscript(commandPart);
      setConversationMode(true);
      sendVoiceCommand(commandPart, true, false);
    } else {
      // Just "Hey Monday" - activate conversation mode
      setConversationMode(true);
      setIsInConversation(true);
      speak("Yes? I'm listening.");
    }
    
    // Reset transcript immediately for next input
    resetTranscript();
  }, [isSpeaking, stopSpeaking, speak, resetTranscript]);

  const handleContinuousConversation = useCallback((fullTranscript: string) => {
    console.log("Processing continuous conversation:", fullTranscript);
    
    if (fullTranscript.trim().length === 0 || isProcessingVoice) {
      return;
    }

    // Stop any ongoing speech immediately
    if (isSpeaking) {
      stopSpeaking();
    }

    // Check for exit commands
    const transcriptLower = fullTranscript.toLowerCase();
    if (transcriptLower.includes("goodbye") || transcriptLower.includes("bye monday") || 
        transcriptLower.includes("stop listening") || transcriptLower.includes("end conversation")) {
      setConversationMode(false);
      setIsInConversation(false);
      speak("Goodbye! Say Hey Monday anytime you want to chat again.");
      resetTranscript();
      return;
    }

    // Process command immediately
    setUserTranscript(fullTranscript);
    sendVoiceCommand(fullTranscript, false, false);
    
    // Reset transcript immediately for next input
    resetTranscript();
  }, [isSpeaking, stopSpeaking, speak, resetTranscript, isProcessingVoice]);

  // Handle incoming WebSocket messages
  useEffect(() => {
    if (wsLastMessage) {
      try {
        const message: WebSocketMessage = JSON.parse(wsLastMessage.data as string);
        console.log("Received WS message:", message);
        
        switch (message.type) {
          case "voice_response":
            handleVoiceResponse(message);
            break;
          case "reasoning_progress":
          case "research_progress":
            handleProgressUpdate(message);
            break;
          case "voice_error":
            setApiError(message.error || "An error occurred");
            setIsThinking(false);
            setIsReceiving(false);
            setIsProcessingVoice(false);
            break;
          case "spatial_response":
            console.log("Spatial response received:", message);
            break;
        }
      } catch (e) {
        console.error("Failed to parse WebSocket message:", e);
        setApiError("Received an invalid message from the server.");
      }
    }
  }, [wsLastMessage]);

  const handleVoiceResponse = (message: WebSocketMessage) => {
    console.log("Handling voice response:", message);
    
    // Immediately clear thinking state when response arrives
    setIsThinking(false);
    setProgressUpdates([]);
    
    // Set the response text for progressive display
    if (message.message) {
      setFullResponseText(message.message);
      setIsReceiving(true);
      setCharIndex(0);
      
      // Speak the response using ElevenLabs TTS
      if (!isSpeaking) {
        speak(message.message);
      }
    }
    
    // Update mode and other data
    if (message.data) {
      if (message.data.mode) {
        setCurrentMode(message.data.mode as "basic" | "reasoning" | "deep-research");
      }
      
      // Handle spatial panels
      if (message.data.panels && message.data.panels.length > 0) {
        setSpatialPanels(message.data.panels);
        
        // Extract visualization data for 3D component
        const mainPanel = message.data.panels.find(p => p.type === 'content');
        if (mainPanel) {
          setVisualizationData({
            query: message.data.query,
            mode: message.data.mode,
            model: message.data.model,
            content: mainPanel.fullContent || mainPanel.content,
            citations: mainPanel.citations || [],
            reasoning: mainPanel.reasoning || []
          });
        }
      }
      
      // Extract YouTube video ID if present
      if (message.data.metadata?.youtubeVideoId) {
        setYoutubeVideoId(message.data.metadata.youtubeVideoId);
      }

      // Check if this was just an activation response
      if (message.data.metadata?.isActivation) {
        // Don't set full conversation state for activation responses
        return;
      }
    }
    
    setIsInConversation(true);
    setIsProcessingVoice(false);
  };

  const handleProgressUpdate = (message: WebSocketMessage) => {
    console.log("Progress update:", message);
    if (message.update) {
      // Only show brief progress updates, don't accumulate them
      setProgressUpdates([message.update]);
      
      // Don't update the current response with progress - keep it separate
      // Progress updates are just for status, not for the final response content
    }
  };

  // Progressive typewriter effect for the main response (significantly sped up)
  useEffect(() => {
    if (isReceiving && fullResponseText && charIndex < fullResponseText.length) {
      const timeout = setTimeout(() => {
        setCurrentResponse(fullResponseText.slice(0, charIndex + 1));
        setCharIndex(charIndex + 1);
      }, 5 + Math.random() * 5); // Reduced from 25-40ms to 5-10ms per character
      return () => clearTimeout(timeout);
    } else if (fullResponseText && charIndex >= fullResponseText.length) {
      setIsReceiving(false);
    }
  }, [charIndex, fullResponseText, isReceiving]);

  // Enhanced voice trigger detection for natural language
  const detectVoiceMode = useCallback((command: string): { mode: "basic" | "reasoning" | "deep-research", confidence: number } => {
    const commandLower = command.toLowerCase().trim();
    
    // Reasoning mode triggers - comprehensive natural language patterns
    const reasoningTriggers = [
      "think about", "analyze", "reasoning", "figure out", "work through", "solve this",
      "explain why", "help me understand", "break down", "think through", "reason about",
      "logic behind", "walk me through", "make sense of", "explain the reasoning",
      "process this", "work this out", "think it through", "analyze this", "help me think",
      "what's the logic", "how does this work", "can you reason", "step by step",
      "logical analysis", "critical thinking", "problem solving", "analytical approach"
    ];
    
    // Deep research mode triggers - comprehensive investigation patterns
    const researchTriggers = [
      "research into", "investigate", "deep dive", "find information about", "look into",
      "study", "explore", "deep research", "comprehensive analysis", "thorough investigation",
      "research this", "dig deeper", "find out about", "learn everything about",
      "comprehensive study", "detailed research", "in-depth analysis", "full investigation",
      "tell me everything", "complete overview", "detailed explanation", "extensive research",
      "academic research", "scholarly analysis", "comprehensive review", "thorough study"
    ];
    
    // Basic mode triggers - quick search and simple queries
    const basicTriggers = [
      "search the web", "find online", "what is", "define", "quick question",
      "search for", "look up", "find", "basic search", "simple query",
      "quick search", "web search", "google", "find me", "search", "who is",
      "when did", "where is", "how much", "how many", "simple answer",
      "quick answer", "basic question", "fast search", "just tell me"
    ];
    
    // Check reasoning triggers first (highest priority for analytical thinking)
    for (const trigger of reasoningTriggers) {
      if (commandLower.includes(trigger)) {
        return { mode: "reasoning", confidence: 0.95 };
      }
    }
    
    // Check research triggers (high priority for comprehensive investigation)
    for (const trigger of researchTriggers) {
      if (commandLower.includes(trigger)) {
        return { mode: "deep-research", confidence: 0.9 };
      }
    }
    
    // Check basic triggers (standard search queries)
    for (const trigger of basicTriggers) {
      if (commandLower.includes(trigger)) {
        return { mode: "basic", confidence: 0.8 };
      }
    }
    
    // Advanced pattern matching for implicit reasoning requests
    if (commandLower.match(/\b(why|how|because|reason|cause|effect|impact|consequence)\b/)) {
      return { mode: "reasoning", confidence: 0.7 };
    }
    
    // Advanced pattern matching for research requests
    if (commandLower.match(/\b(history|background|origin|development|evolution|comprehensive|detailed)\b/)) {
      return { mode: "deep-research", confidence: 0.6 };
    }
    
    // Default to current mode with low confidence
    return { mode: currentMode, confidence: 0.1 };
  }, [currentMode]);

  // Enhanced mode switching with visual feedback
  const switchToMode = useCallback((newMode: "basic" | "reasoning" | "deep-research", command: string, isVoiceTriggered: boolean = true) => {
    if (newMode === currentMode) return;
    
    console.log(`Voice-triggered mode switch: ${currentMode} → ${newMode}`);
    
    // Set immediate visual feedback
    setIsModeSwitching(true);
    
    // Create mode-specific status message
    const modeLabels = {
      "basic": "Basic Search",
      "reasoning": "Reasoning Mode", 
      "deep-research": "Deep Research Mode"
    };
    
    const switchMessage = `Switching to ${modeLabels[newMode]}...`;
    setModeSwitchMessage(switchMessage);
    
    // Immediate mode change for instant visual feedback
    setCurrentMode(newMode);
    
    // Voice confirmation of mode switch
    if (isVoiceTriggered && !isSpeaking) {
      speak(`Switching to ${modeLabels[newMode]}`);
    }
    
    // Clear switch message after brief display
    setTimeout(() => {
      setModeSwitchMessage("");
      setIsModeSwitching(false);
    }, 2000);
    
  }, [currentMode, isSpeaking, speak]);

  const sendVoiceCommand = useCallback((command: string, isExplicitTrigger: boolean = false, isActivation: boolean = false) => {
    if (!wsIsConnected) {
      setApiError("Not connected to Monday backend. Please wait for connection.");
      return;
    }

    if (isProcessingVoice && !isExplicitTrigger) {
      console.log("Already processing voice, ignoring command");
      return;
    }

    console.log("Sending voice command:", command);
    
    // FIRST: Detect and switch mode immediately for instant visual feedback
    const modeDetection = detectVoiceMode(command);
    if (modeDetection.confidence > 0.5) {
      switchToMode(modeDetection.mode, command, true);
    }
    
    // Clear previous state immediately
    setApiError(null);
    setIsThinking(true);
    setIsReceiving(false);
    setIsProcessingVoice(true);
    setCurrentResponse("");
    setFullResponseText("");
    setCharIndex(0);
    setProgressUpdates([]);

    // Send command via WebSocket immediately
    wsSendMessage(JSON.stringify({
      type: "voice_command",
      command: command,
      conversationActive: isInConversation,
      isExplicitTrigger: isExplicitTrigger,
      isActivation: isActivation,
      timestamp: Date.now()
    }));
  }, [wsIsConnected, wsSendMessage, isInConversation, isProcessingVoice, detectVoiceMode, switchToMode]);

  const toggleListening = () => {
    if (!isMounted) return;
    
    if (listening) {
      SpeechRecognition.stopListening();
      setIsListeningForTrigger(false);
      setConversationMode(false);
    } else {
      startAutomaticListening();
    }
  };

  const startNewResponseFlow = (mode: "basic" | "reasoning" | "deep-research") => {
    setCurrentMode(mode);
    setUserTranscript("");
    setCurrentResponse("");
    setFullResponseText("");
    setCharIndex(0);
    setIsThinking(false);
    setIsReceiving(false);
    setApiError(null);
    setSpatialPanels([]);
    setProgressUpdates([]);
    setIsProcessingVoice(false);
    setConversationMode(false);
    
    // Stop any ongoing speech
    if (isSpeaking) {
      stopSpeaking();
    }
  };

  const getModeConfig = (mode: string) => {
    const configs = {
      "reasoning": { 
        icon: <Brain className="w-4 h-4" />, 
        label: "Reasoning Pro",
        color: "from-purple-500 to-purple-700",
        description: "Deep analytical thinking"
      },
      "deep-research": { 
        icon: <Activity className="w-4 h-4" />, 
        label: "Deep Research",
        color: "from-blue-500 to-blue-700", 
        description: "Comprehensive investigation"
      },
      "basic": { 
        icon: <Search className="w-4 h-4" />, 
        label: "Basic",
        color: "from-green-500 to-green-700",
        description: "Quick search & answers"
      }
    };
    
    return configs[mode as keyof typeof configs] || configs.basic;
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] p-8 text-[#FBFAF4]">
      <div className="w-full max-w-[2000px] mx-auto">
        <div className="text-center mb-12">
          <div className="flex items-center justify-center gap-6 mb-6">
            <PerplexityLogo className="w-16 h-16" />
            <div className="text-center">
              <h1 className="text-5xl font-bold tracking-tight mb-2">Monday</h1>
              <p className="text-[#20808D] text-lg font-semibold">Learning, Amalgamated</p>
            </div>
            <div className="flex items-center gap-3 bg-gradient-to-r from-[#20808D]/20 to-[#20808D]/10 rounded-2xl px-6 py-3 backdrop-blur-xl border border-[#20808D]/30">
              <Zap className="w-5 h-5 text-[#20808D]" />
              <span className="text-[#20808D] text-sm font-bold">VR AI Companion</span>
              <Sparkles className="w-4 h-4 text-[#20808D] animate-pulse" />
            </div>
          </div>
          
          {/* Connection Status */}
          <div className="flex items-center justify-center gap-2 mb-4">
            <div className={`w-3 h-3 rounded-full ${wsIsConnected ? 'bg-green-400 animate-pulse' : 'bg-red-400'}`}></div>
            <span className="text-sm text-[#20808D]">
              {wsIsConnected ? 'Connected to Monday Backend' : 'Connecting to Monday Backend...'}
            </span>
          </div>
        </div>

        {/* Enhanced Mode Switching UI with Visual Feedback */}
        <div className="flex justify-center gap-4 mb-12">
          {(["basic", "reasoning", "deep-research"] as const).map((mode) => {
            const config = getModeConfig(mode);
            const isActive = currentMode === mode;
            return (
              <button
                key={mode}
                onClick={() => startNewResponseFlow(mode)}
                className={`relative flex items-center gap-3 px-6 py-3 rounded-2xl transition-all duration-300 text-sm font-bold backdrop-blur-xl border transform hover:scale-105 ${
                  isActive
                    ? `bg-gradient-to-r ${config.color} text-white border-white/50 shadow-xl shadow-current/30 scale-110`
                    : "bg-gradient-to-r from-[#20808D]/20 to-[#20808D]/10 text-[#20808D] hover:from-[#20808D]/30 hover:to-[#20808D]/20 border-[#20808D]/30"
                }`}
              >
                {isActive && isModeSwitching && (
                  <div className="absolute inset-0 rounded-2xl bg-gradient-to-r from-white/20 to-transparent animate-pulse"></div>
                )}
                {config.icon}
                <div className="flex flex-col items-start">
                  <span>Sonar {config.label}</span>
                  {isActive && (
                    <span className="text-xs opacity-80">{config.description}</span>
                  )}
                </div>
                {isActive && (
                  <div className="absolute -top-1 -right-1 w-3 h-3 bg-white rounded-full animate-pulse"></div>
                )}
              </button>
            );
          })}
        </div>

        {/* Mode Switch Status Message */}
        {modeSwitchMessage && (
          <div className="text-center mb-4 p-3 bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-500/30 rounded-xl">
            <div className="flex items-center justify-center gap-2">
              <div className="w-4 h-4 border-2 border-blue-400 border-t-transparent rounded-full animate-spin"></div>
              <span className="text-blue-400 font-semibold">{modeSwitchMessage}</span>
            </div>
          </div>
        )}

        {/* Enhanced Status Indicators with Mode Information */}
        {(isSpeaking || isProcessingVoice || isModeSwitching) && (
          <div className="text-center mb-4 text-[#20808D] text-sm">
            {isSpeaking && <span className="mr-4">🔊 Monday is speaking...</span>}
            {isProcessingVoice && (
              <span className="mr-4">
                {currentMode === 'reasoning' && '🧠 Analyzing with deep reasoning...'}
                {currentMode === 'deep-research' && '🔍 Conducting comprehensive research...'}
                {currentMode === 'basic' && '🤔 Processing your request...'}
              </span>
            )}
            {isModeSwitching && <span>⚡ Switching modes...</span>}
          </div>
        )}

        {/* Enhanced Voice Status with Mode Display */} 
        <div className="flex justify-center items-center gap-4 mb-8">
          <button 
            onClick={toggleListening}
            disabled={!isMounted || !browserSupportsSpeechRecognition || !isMicrophoneAvailable}
            className={`p-4 rounded-full transition-colors duration-200 ease-in-out shadow-xl group focus:outline-none focus:ring-2 focus:ring-[#20808D] focus:ring-opacity-50 ${
              listening ? 'bg-green-500 hover:bg-green-600' : 'bg-gray-500 hover:bg-gray-600'
            } ${(!isMounted || !browserSupportsSpeechRecognition || !isMicrophoneAvailable) ? 'opacity-50 cursor-not-allowed' : ''}`}
            title={listening ? "Voice Recognition Active" : "Voice Recognition Inactive"}
          >
            {listening ? <Mic className="w-8 h-8 text-white animate-pulse" /> : <MicOff className="w-8 h-8 text-white" />}
          </button>
          
          <div className="text-center">
            <div className="text-sm text-[#20808D] mb-1">
              {conversationMode 
                ? `💬 In conversation - ${getModeConfig(currentMode).label} mode active`
                : listening 
                ? "🎤 Listening for 'Hey Monday'..." 
                : "Voice recognition inactive"
              }
            </div>
            {!conversationMode && isListeningForTrigger && (
              <div className="text-xs text-gray-400">
                Say "Hey Monday" + command to start (auto-detects mode)
              </div>
            )}
            {conversationMode && (
              <div className="text-xs text-gray-400">
                Voice triggers: "think about", "research into", "search for"
              </div>
            )}
          </div>
        </div>
        
        {/* Error Display */}
        {apiError && (
          <div className="text-center text-red-400 mb-4 p-4 bg-red-500/10 border border-red-500/30 rounded-xl">
            {apiError}
          </div>
        )}
        
        {/* User Transcript Display */}
        {userTranscript && (
          <div className="text-center text-gray-400 mb-4 p-4 bg-[#20808D]/10 border border-[#20808D]/20 rounded-xl">
            <span className="font-semibold text-gray-300">You said:</span> {userTranscript}
          </div>
        )}

        <div className="grid grid-cols-1 xl:grid-cols-12 gap-8 h-[800px]">
          <div className="xl:col-span-4 transform xl:-rotate-2 xl:scale-95 origin-right">
            <div className="h-full shadow-2xl shadow-[#20808D]/20">
              <BinaryTree3D visualizationData={visualizationData} title="Query Visualization" /> 
            </div>
          </div>

          <div className="xl:col-span-4 z-10">
            <div className="h-full shadow-2xl shadow-[#20808D]/30">
              <VoiceProcessingPanel
                currentResponse={currentResponse}
                isThinking={isThinking}
                isTyping={isReceiving}
                mode={currentMode}
                userTranscript={userTranscript}
                voiceActive={listening}
                processingError={apiError}
                isModeSwitching={isModeSwitching}
                modeSwitchMessage={modeSwitchMessage}
              />
            </div>
          </div>

          <div className="xl:col-span-4 transform xl:rotate-2 xl:scale-95 origin-left">
            <div className="h-full shadow-2xl shadow-[#20808D]/20">
              <YouTubePanel 
                videoId={youtubeVideoId} 
                query={visualizationData?.query || userTranscript}
                mode={currentMode}
              />
            </div>
          </div>
        </div>

        <div className="text-center mt-12">
          <div className="inline-flex items-center gap-4 bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 rounded-2xl px-8 py-4 backdrop-blur-xl border border-[#20808D]/20">
            <PerplexityLogo className="w-6 h-6" />
            <span className="text-[#20808D]/80 text-sm font-medium">
              Powered by Perplexity Sonar API • ElevenLabs TTS • WebXR Ready • Quest 2/3 Compatible
            </span>
          </div>
        </div>
      </div>
    </div>
  );
}
