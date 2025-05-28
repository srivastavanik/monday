"use client"

import { useState, useEffect, useCallback, useRef } from "react"
import { AdaptiveVisualizationPanel } from "./components/adaptive-visualization-panel"
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
  type: "voice_response" | "reasoning_progress" | "research_progress" | "voice_error" | "spatial_response" | "heartbeat_ack";
  message?: string;
  data?: {
    panels?: SpatialPanel[];
    mode?: string;
    model?: string;
    query?: string;
    citations?: any[];
    reasoning?: any[];
    sources?: any[];
    thinking?: string;
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
  // --- Group 1: State Hooks ---
  const [currentMode, setCurrentMode] = useState<"basic" | "reasoning" | "deep-research">("basic");
  const [currentResponse, setCurrentResponse] = useState("");
  const [fullResponseText, setFullResponseText] = useState("");
  const [isThinking, setIsThinking] = useState(false);
  const [isReceiving, setIsReceiving] = useState(false);
  const [charIndex, setCharIndex] = useState(0);
  const [youtubeVideoId, setYoutubeVideoId] = useState<string | undefined>(undefined);
  const [visualizationData, setVisualizationData] = useState<any>(null);
  const [userTranscript, setUserTranscript] = useState("");
  const [apiError, setApiError] = useState<string | null>(null);
  const [isMounted, setIsMounted] = useState(false);
  const [spatialPanels, setSpatialPanels] = useState<SpatialPanel[]>([]);
  const [isInConversation, setIsInConversation] = useState(false);
  const [progressUpdates, setProgressUpdates] = useState<string[]>([]);
  const [isProcessingVoice, setIsProcessingVoice] = useState(false);
  const [isListeningForTrigger, setIsListeningForTrigger] = useState(true);
  const [lastTriggerTime, setLastTriggerTime] = useState(0);
  const [conversationMode, setConversationMode] = useState(false);
  const [modeSwitchMessage, setModeSwitchMessage] = useState<string>("");
  const [isModeSwitching, setIsModeSwitching] = useState(false);
  const [liveTranscript, setLiveTranscript] = useState<string>("");
  const [lastProcessedTranscript, setLastProcessedTranscript] = useState<string>("");
  const [progressSources, setProgressSources] = useState<any[]>([]);
  const [progressReasoning, setProgressReasoning] = useState<any[]>([]);
  const [progressCitations, setProgressCitations] = useState<any[]>([]);

  // --- Group 2: Custom Hooks ---
  const { speak, speakProgress, stop: stopSpeaking, isSpeaking, initializeAudioContext, isInitialized } = useTextToSpeech({
    voiceId: "EXAVITQu4vr4xnSDxMaL", // Sarah - available in your account
    apiKey: "sk_6be41ddad8e1f62baabff7344e221588d655f06c35823991"
  });

  const {
    transcript,
    listening,
    resetTranscript,
    browserSupportsSpeechRecognition,
    isMicrophoneAvailable
  } = useSpeechRecognition();

  // --- Group 3: WebSocket Hook Setup ---
  const handleWebSocketOpen = useCallback(() => {
    console.log("FRONTEND: WebSocket connected successfully.");
    setApiError(null); 
  }, []);

  const handleWebSocketClose = useCallback((event: CloseEvent) => {
    console.warn("FRONTEND: WebSocket disconnected.", "Reason:", event.reason, "Code:", event.code, "WasClean:", event.wasClean);
    if (isMounted) { 
        if (!event.wasClean) {
            setApiError("WebSocket connection lost. Attempting to reconnect...");
        } else {
            setApiError("WebSocket connection closed."); 
        }
    }
  }, [isMounted]); 

  const handleWebSocketError = useCallback((event: Event) => {
    console.error("FRONTEND: WebSocket connection error event:", event);
    if (isMounted) { 
        setApiError("WebSocket connection error. Check backend and network.");
    }
  }, [isMounted]); 

  const { 
    isConnected: wsIsConnected,
    lastMessage: wsLastMessage,
    sendMessage: wsSendMessage,
    error: wsHookError,
    reconnect: wsReconnect 
  } = useMondayWebSocket({
    socketUrl: MONDAY_WEBSOCKET_URL,
    onOpen: handleWebSocketOpen,
    onClose: handleWebSocketClose, 
    onError: handleWebSocketError,
  });

  // --- Group 4: Core Voice Processing Functions (in dependency order) ---
  const detectVoiceMode = useCallback((command: string): { mode: "basic" | "reasoning" | "deep-research", confidence: number } => {
    const commandLower = command.toLowerCase().trim();
    const reasoningTriggers = ["think about", "analyze", "reasoning", "figure out", "work through", "solve this", "explain why", "help me understand", "break down", "think through", "reason about", "logic behind", "walk me through", "make sense of", "explain the reasoning", "process this", "work this out", "think it through", "analyze this", "help me think", "what's the logic", "how does this work", "can you reason", "step by step", "logical analysis", "critical thinking", "problem solving", "analytical approach"];
    const researchTriggers = ["research into", "investigate", "deep dive", "find information about", "look into", "study", "explore", "deep research", "comprehensive analysis", "thorough investigation", "research this", "dig deeper", "find out about", "learn everything about", "comprehensive study", "detailed research", "in-depth analysis", "full investigation", "tell me everything", "complete overview", "detailed explanation", "extensive research", "academic research", "scholarly analysis", "comprehensive review", "thorough study"];
    const basicTriggers = ["search the web", "find online", "what is", "define", "quick question", "search for", "look up", "find", "basic search", "simple query", "quick search", "web search", "google", "find me", "search", "who is", "when did", "where is", "how much", "how many", "simple answer", "quick answer", "basic question", "fast search", "just tell me"];
    for (const trigger of reasoningTriggers) if (commandLower.includes(trigger)) return { mode: "reasoning", confidence: 0.95 };
    for (const trigger of researchTriggers) if (commandLower.includes(trigger)) return { mode: "deep-research", confidence: 0.9 };
    for (const trigger of basicTriggers) if (commandLower.includes(trigger)) return { mode: "basic", confidence: 0.8 };
    if (commandLower.match(/\b(why|how|because|reason|cause|effect|impact|consequence)\b/)) return { mode: "reasoning", confidence: 0.7 };
    if (commandLower.match(/\b(history|background|origin|development|evolution|comprehensive|detailed)\b/)) return { mode: "deep-research", confidence: 0.6 };
    return { mode: currentMode, confidence: 0.1 };
  }, [currentMode]);

  // Handle YouTube playback commands (moved before useEffect)
  const handleYouTubeCommand = useCallback((command: string) => {
    const youtubeFrame = document.getElementById('youtube-player');
    if (!youtubeFrame) {
      console.log("YouTube player not found");
      return;
    }

    try {
      const frame = youtubeFrame as HTMLIFrameElement;
      if (command.includes('play')) {
        frame.contentWindow?.postMessage('{"event":"command","func":"playVideo","args":""}', '*');
        speak("Playing video");
      } else if (command.includes('pause') || command.includes('stop')) {
        frame.contentWindow?.postMessage('{"event":"command","func":"pauseVideo","args":""}', '*');
        speak("Pausing video");
      } else if (command.includes('mute')) {
        frame.contentWindow?.postMessage('{"event":"command","func":"mute","args":""}', '*');
        speak("Muting video");
      } else if (command.includes('unmute')) {
        frame.contentWindow?.postMessage('{"event":"command","func":"unMute","args":""}', '*');
        speak("Unmuting video");
      } else if (command.includes('restart')) {
        frame.contentWindow?.postMessage('{"event":"command","func":"seekTo","args":[0, true]}', '*');
        speak("Restarting video");
      }
    } catch (error) {
      console.error("Error controlling YouTube player:", error);
    }
  }, [speak]);

  const switchToMode = useCallback((newMode: "basic" | "reasoning" | "deep-research", command: string, isVoiceTriggered: boolean = true) => {
    if (newMode === currentMode) return;
    setIsModeSwitching(true);
    const modeLabels = { "basic": "Basic Search", "reasoning": "Reasoning Mode", "deep-research": "Deep Research Mode" };
    setModeSwitchMessage(`Switching to ${modeLabels[newMode]}...`);
    setCurrentMode(newMode);
    if (isVoiceTriggered && !isSpeaking) speak(`Switching to ${modeLabels[newMode]}`);
    setTimeout(() => { setModeSwitchMessage(""); setIsModeSwitching(false); }, 2000);
  }, [currentMode, isSpeaking, speak]);

  const sendVoiceCommand = useCallback((command: string, isExplicitTrigger: boolean = false, isActivation: boolean = false) => {
    if (!command.trim() || (isProcessingVoice && !isExplicitTrigger)) {
      console.log("Skipping voice command - empty or already processing:", { command, isProcessingVoice, isExplicitTrigger });
      return;
    }
    
    if (!wsIsConnected) {
      console.error("Cannot send voice command - WebSocket not connected");
      setApiError("Cannot send command: Backend not connected. Please try again.");
      return;
    }

    const modeDetection = detectVoiceMode(command);
    let actualMode = currentMode;
    if (modeDetection.confidence > 0.5) {
      switchToMode(modeDetection.mode, command, true);
      actualMode = modeDetection.mode;
    }
    
    setApiError(null);
    // Only set thinking state for reasoning and deep-research modes
    if (actualMode === 'reasoning' || actualMode === 'deep-research') {
      setIsThinking(true);
      setIsProcessingVoice(true);
    } else {
      // For basic mode, don't show processing state
      setIsThinking(false);
      setIsProcessingVoice(false);
      // Immediately set visualization data for basic mode
      setVisualizationData({
        query: command,
        mode: actualMode,
        model: 'sonar',
        content: '',
        citations: []
      });
    }
    
    setIsReceiving(false); 
    setCurrentResponse(""); 
    setFullResponseText(""); 
    setCharIndex(0);
    setProgressUpdates([]); 
    setProgressSources([]); 
    setProgressReasoning([]); 
    setProgressCitations([]);
    setLiveTranscript("");
    
    console.log("Sending voice command:", { command, mode: actualMode, isExplicitTrigger, isActivation });
    
    try {
      const commandData = {
        type: "voice_command", 
        command: command, 
        conversationActive: isInConversation, 
        isExplicitTrigger: isExplicitTrigger, 
        isActivation: isActivation, 
        timestamp: Date.now()
      };
      
      wsSendMessage(JSON.stringify(commandData));
      console.log("Voice command sent successfully");
    } catch (error) {
      console.error("Failed to send voice command:", error);
      setApiError("Failed to send command. Please check connection.");
      setIsThinking(false); 
      setIsProcessingVoice(false);
    }
  }, [wsSendMessage, isInConversation, isProcessingVoice, detectVoiceMode, switchToMode, currentMode, wsIsConnected]);

  const startAutomaticListening = useCallback(() => {
    if (!browserSupportsSpeechRecognition) {
      setApiError("Voice recognition not supported in this browser."); 
      return;
    }
    
    if (!isMicrophoneAvailable) {
      console.log("Microphone not available, requesting permissions...");
      // Try to get permissions explicitly
      navigator.mediaDevices.getUserMedia({ audio: true })
        .then(stream => {
          console.log("Microphone permission granted, stopping test stream...");
          stream.getTracks().forEach(track => track.stop());
          // Try again after getting permission
          setTimeout(() => {
            if (isMounted && !conversationMode && !listening) {
              startAutomaticListening();
            }
          }, 1000);
        })
        .catch(err => {
          console.error("Microphone permission denied:", err);
          setApiError("Microphone access required for hands-free operation. Please allow permissions and refresh.");
        });
      return;
    }
    
    if (conversationMode || listening) {
      console.log("Already in conversation mode or listening, skipping...");
      return;
    }
    
    console.log("🎤 Starting hands-free speech recognition...");
    console.log("Current state:", { 
      isMounted, 
      browserSupportsSpeechRecognition, 
      isMicrophoneAvailable, 
      listening, 
      conversationMode 
    });
    
    try {
      // Use abort controller for better control
      const abortController = new AbortController();
      
      SpeechRecognition.startListening({ 
        continuous: true, 
        language: 'en-US', 
        interimResults: true 
      }).then(() => {
        console.log("✅ Hands-free voice recognition started successfully!");
        setIsListeningForTrigger(true);
        setApiError(null);
      }).catch((error: any) => {
        console.error("❌ Speech recognition failed to start:", error);
        
        // Retry with a delay if it's a temporary issue
        if (error.message && error.message.includes("already started")) {
          console.log("Speech recognition already running, that's fine");
          setIsListeningForTrigger(true);
        } else {
          console.log("Will retry hands-free activation in 2 seconds...");
          setTimeout(() => {
            if (isMounted && !conversationMode && !listening) {
              startAutomaticListening();
            }
          }, 2000);
        }
      });
      
    } catch (error) {
      console.error("Error in speech recognition setup:", error);
      // Retry after a delay
      setTimeout(() => {
        if (isMounted && !conversationMode && !listening) {
          console.log("Retrying hands-free voice activation...");
          startAutomaticListening();
        }
      }, 2000);
    }
  }, [browserSupportsSpeechRecognition, isMicrophoneAvailable, listening, conversationMode, isMounted]);

  const handleMondayActivation = useCallback((fullTranscript: string) => {
    if (isSpeaking) stopSpeaking();
    const triggerPattern = /(hey monday|hey monty)/gi;
    const matches = [...fullTranscript.matchAll(triggerPattern)];
    if (matches.length === 0) return;
    const lastMatch = matches[matches.length - 1];
    const commandPart = fullTranscript.substring((lastMatch.index || 0) + lastMatch[0].length).trim();
    const cleanCommand = commandPart.replace(triggerPattern, '').trim();
    
    setConversationMode(true); 
    setIsInConversation(true); 
    setLiveTranscript(""); 
    resetTranscript();
    setIsListeningForTrigger(false);

    if (cleanCommand) {
      setUserTranscript(cleanCommand);
      sendVoiceCommand(cleanCommand, true, false);
    } else {
      sendVoiceCommand("", true, true);
    }
  }, [isSpeaking, stopSpeaking, resetTranscript, sendVoiceCommand]);

  const handleContinuousConversation = useCallback((fullTranscript: string) => {
    if (fullTranscript.trim().length === 0 || isProcessingVoice) return;
    
    const transcriptLower = fullTranscript.toLowerCase();
    if (transcriptLower.includes("goodbye") || transcriptLower.includes("bye monday") || transcriptLower.includes("stop listening") || transcriptLower.includes("end conversation")) {
      setConversationMode(false); 
      setIsInConversation(false); 
      setLiveTranscript(""); 
      resetTranscript();
      speak("Goodbye! Say Hey Monday anytime you want to chat again.");
      return;
    }
    
    const cleanTranscript = fullTranscript.replace(/(hey monday|hey monty)/gi, '').replace(/\s+/g, ' ').trim();
    if (cleanTranscript.length > 0 && cleanTranscript !== lastProcessedTranscript) {
      console.log("Processing voice command:", cleanTranscript);
      setLastProcessedTranscript(cleanTranscript);
      setUserTranscript(cleanTranscript);
      
      setTimeout(() => {
        if (wsIsConnected) {
          sendVoiceCommand(cleanTranscript, false, false);
        } else {
          setApiError("Cannot send command: Backend connection lost");
          console.error("WebSocket not connected, cannot send voice command");
        }
      }, 100);
    }
  }, [isProcessingVoice, resetTranscript, lastProcessedTranscript, speak, sendVoiceCommand, wsIsConnected]);

  const handleVoiceResponse = useCallback((message: WebSocketMessage) => {
    console.log("HandleVoiceResponse called with:", { 
      type: message.type, 
      hasMessage: !!message.message,
      messagePreview: message.message?.substring(0, 50),
      isActivation: message.data?.metadata?.isActivation,
      mode: message.data?.mode
    });
    
    if (message.data?.metadata?.isActivation && message.message) {
      console.log("Handling activation greeting:", message.message);
      setFullResponseText(message.message);
      setIsReceiving(true);
      setCharIndex(0);
      setIsProcessingVoice(false);
      
      setTimeout(() => { 
        if (!isSpeaking && message.message) {
          console.log("Speaking activation greeting...");
          speak(message.message);
        }
      }, 500);
      return;
    }
    
    // For basic mode, immediately clear processing states
    if (message.data?.mode === 'basic') {
      setIsThinking(false);
      setIsProcessingVoice(false);
      setIsReceiving(true);
    }
    
    if (message.message) {
      setFullResponseText(message.message);
      setIsReceiving(true);
      setCharIndex(0);
      
      // Clear thinking state for all modes when response arrives
      if (message.data && (message.data.mode === 'reasoning' || message.data.mode === 'deep-research')) {
        setIsThinking(false);
        setIsProcessingVoice(false);
      }
      
      setTimeout(() => { 
        if (!isSpeaking && message.message) {
          console.log("Speaking response...");
          speak(message.message);
        }
      }, 500);
    }
    
    if (message.data) {
      if (message.data.mode) setCurrentMode(message.data.mode as "basic" | "reasoning" | "deep-research");
      if (message.data.query && message.data.query !== 'activation') setUserTranscript(message.data.query);
      if (message.data.sources?.length) setProgressSources(message.data.sources);
      if (message.data.citations?.length) setProgressCitations(message.data.citations);
      
      // Update visualization data for all modes
      if (message.data.query && message.data.query !== 'activation') {
        setVisualizationData({ 
          query: message.data.query, 
          mode: message.data.mode, 
          model: message.data.model || 'sonar',
          content: message.message || message.data.query,
          citations: message.data.citations || [],
          reasoning: message.data.reasoning || [],
          thinking: message.data.thinking || ''
        });
      }
      
      if (message.data.panels?.length) {
        setSpatialPanels(message.data.panels);
        const mainPanel = message.data.panels.find(p => p.type === 'content');
        if (mainPanel) {
          setVisualizationData({ 
            query: message.data.query, 
            mode: message.data.mode, 
            model: message.data.model, 
            content: mainPanel.fullContent || mainPanel.content, 
            citations: mainPanel.citations || [], 
            reasoning: mainPanel.reasoning || [],
            thinking: message.data.thinking || ''
          });
        }
      }
      
      if (message.data.metadata?.youtubeVideoId) {
        console.log("Setting YouTube video ID:", message.data.metadata.youtubeVideoId);
        setYoutubeVideoId(message.data.metadata.youtubeVideoId);
      } else if (message.data.mode === 'basic' && message.data.query && message.data.query !== 'activation') {
        console.log("Resetting YouTube video ID for new query");
        setYoutubeVideoId(undefined);
      }
    }
    
    setIsInConversation(true);
    setIsProcessingVoice(false);
  }, [isSpeaking, speak]);

  const handleProgressUpdate = useCallback((message: WebSocketMessage) => {
    if (message.update) {
      setProgressUpdates(prev => [...prev, message.update as string]);
      if (message.data) {
        if (message.data.sources?.length) setProgressSources(message.data.sources);
        if (message.data.reasoning?.length) setProgressReasoning(message.data.reasoning);
        if (message.data.citations?.length) setProgressCitations(message.data.citations);
      }
    }
  }, []);

  const toggleListening = useCallback(() => {
    if (!isMounted) return;
    if (listening || conversationMode) {
      SpeechRecognition.stopListening();
      setIsListeningForTrigger(false); 
      setConversationMode(false); 
      setIsInConversation(false);
      setLiveTranscript(""); 
      resetTranscript();
    } else {
      if (isMicrophoneAvailable) {
        startAutomaticListening();
      } else {
        setApiError("Microphone not available. Please grant permission first.");
      }
    }
  }, [isMounted, listening, conversationMode, resetTranscript, startAutomaticListening, isMicrophoneAvailable]);

  const startNewResponseFlow = useCallback((mode: "basic" | "reasoning" | "deep-research") => {
    setCurrentMode(mode); setUserTranscript(""); setCurrentResponse(""); setFullResponseText("");
    setCharIndex(0); setIsThinking(false); setIsReceiving(false); setApiError(null);
    setSpatialPanels([]); setProgressUpdates([]); setProgressSources([]); setProgressReasoning([]);
    setProgressCitations([]); setIsProcessingVoice(false); setConversationMode(false);
    if (isSpeaking) stopSpeaking();
  }, [isSpeaking, stopSpeaking]);

  // --- Group 5: useEffect Hooks ---
  useEffect(() => { 
    setIsMounted(true);
    return () => {
        console.log("MondayPerplexitySystem: Component unmounting. Setting isMounted to false.");
        setIsMounted(false);
    };
  }, []);

  // WebSocket Connection Status/Error Effect
  useEffect(() => {
    if (wsHookError && isMounted) {
        console.error("FRONTEND: WebSocket hook reported an error:", wsHookError);
        let errorMessage = "WebSocket issue. Check backend.";
        if (wsHookError instanceof Error) {
            errorMessage = `WebSocket connection error: ${wsHookError.message}`;
        } else if (typeof wsHookError === 'string') {
            errorMessage = wsHookError;
        } else if (wsHookError.type) {
            errorMessage = `WebSocket event error: ${wsHookError.type}`;
        }
        setApiError(errorMessage);
    }
  }, [wsHookError, isMounted]);

  // Polyfill Effect
  useEffect(() => {
    if (isMounted && !SpeechRecognition.browserSupportsSpeechRecognition() && typeof window !== 'undefined') {
      console.log("Speech recognition not natively supported, attempting polyfill...");
      
      const SpeechRecognitionPolyfill = (window as any).webkitSpeechRecognition || (window as any).SpeechRecognition;
      if (SpeechRecognitionPolyfill) {
        console.log("Found speech recognition polyfill:", SpeechRecognitionPolyfill.name);
        try { 
          SpeechRecognition.applyPolyfill(SpeechRecognitionPolyfill); 
          console.log("Speech recognition polyfill applied successfully");
        } catch (error) { 
          console.error("Polyfill application error:", error); 
          setApiError("Speech recognition initialization failed. Please try refreshing the page.");
        }
      } else {
        console.error("No speech recognition support found in browser");
        setApiError("Speech recognition not supported in this browser. Please use Chrome, Edge, or Safari.");
      }
    } else if (isMounted && SpeechRecognition.browserSupportsSpeechRecognition()) {
      console.log("Native speech recognition support detected");
    }
  }, [isMounted]);

  // Microphone Permission and Audio Context Initialization Effect
  useEffect(() => {
    if (!isMounted || !browserSupportsSpeechRecognition) return;
    initializeAudioContext(); 
    if (isMicrophoneAvailable === false) { 
      console.log("Mic Permission/Audio Init Effect: Microphone NOT available or permission not yet determined. Requesting permissions...");
      navigator.mediaDevices.getUserMedia({ audio: true })
        .then(stream => {
          console.log("Mic Permission/Audio Init Effect: Microphone permission GRANTED by user.");
          stream.getTracks().forEach(track => track.stop());
        })
        .catch(err => {
          console.error("Mic Permission/Audio Init Effect: Microphone permission DENIED by user:", err);
          setApiError("Microphone access denied. Please allow permissions and refresh.");
        });
    } else if (isMicrophoneAvailable === true) {
        console.log("Mic Permission/Audio Init Effect: Microphone is already available.");
    }
  }, [isMounted, browserSupportsSpeechRecognition, isMicrophoneAvailable, initializeAudioContext]);

  // Effect to Start Automatic Listening - RESTORED FOR HANDS-FREE OPERATION
  useEffect(() => {
    const micTrulyAvailable = isMicrophoneAvailable === true;
    console.log("Auto-start voice recognition check:", { 
      isMounted, 
      browserSupportsSpeechRecognition, 
      micTrulyAvailable, 
      listening, 
      conversationMode 
    });
    
    if (isMounted && browserSupportsSpeechRecognition && micTrulyAvailable && !listening && !conversationMode) {
      console.log("Starting hands-free voice recognition automatically...");
      
      // Start immediately for hands-free operation
      const startTimer = setTimeout(() => {
        if (isMounted && !listening && !conversationMode) {
          console.log("Initiating hands-free voice activation...");
          startAutomaticListening();
        }
      }, 500); // Small delay to ensure everything is ready
      
      return () => clearTimeout(startTimer);
    }
  }, [isMounted, browserSupportsSpeechRecognition, isMicrophoneAvailable, listening, conversationMode, startAutomaticListening]);
  
  // Debug Speech State Effect 
  useEffect(() => {
    console.log("Speech State Change:", { 
      browserSupportsSpeechRecognition, 
      isMicrophoneAvailable, 
      listening, 
      transcriptLength: transcript ? transcript.length : 0,
      transcriptPreview: transcript?.substring(0,30),
      conversationMode, 
      isListeningForTrigger,
      wsConnected: wsIsConnected
    });
  }, [browserSupportsSpeechRecognition, isMicrophoneAvailable, listening, transcript, conversationMode, isListeningForTrigger, wsIsConnected]);

  // Improved speech recognition initialization with auto-restart
  useEffect(() => {
    // Handle speech recognition lifecycle for hands-free operation
    if (isMounted && browserSupportsSpeechRecognition && isMicrophoneAvailable === true) {
      console.log('Setting up hands-free speech recognition lifecycle...');
      
      let restartTimer: NodeJS.Timeout;
      
      // Monitor speech recognition state
      const checkAndRestart = () => {
        if (isMounted && !conversationMode && !listening && isListeningForTrigger) {
          console.log('Speech recognition stopped unexpectedly, restarting for hands-free operation...');
          clearTimeout(restartTimer);
          restartTimer = setTimeout(() => {
            if (isMounted && !conversationMode && !listening) {
              startAutomaticListening();
            }
          }, 1500);
        }
      };
      
      // Check periodically if we should be listening but aren't
      const monitorInterval = setInterval(() => {
        if (isMounted && !conversationMode && !listening && isListeningForTrigger) {
          console.log('Voice recognition monitor: Not listening when we should be, restarting...');
          startAutomaticListening();
        }
      }, 5000); // Check every 5 seconds
      
      return () => {
        clearInterval(monitorInterval);
        clearTimeout(restartTimer);
        console.log('Speech recognition lifecycle cleanup');
      };
    }
  }, [isMounted, browserSupportsSpeechRecognition, isMicrophoneAvailable, conversationMode, isListeningForTrigger, listening, startAutomaticListening]);

  // Transcript Processing for "Hey Monday" Trigger Effect
  useEffect(() => {
    if (transcript && listening && !conversationMode && isListeningForTrigger) {
      console.log("=== TRANSCRIPT PROCESSING ===");
      console.log("Raw transcript:", transcript);
      console.log("Transcript length:", transcript.length);
      console.log("Listening state:", listening);
      console.log("Conversation mode:", conversationMode);
      console.log("Is listening for trigger:", isListeningForTrigger);
      
      const transcriptLower = transcript.toLowerCase().trim();
      console.log("Processed transcript (lowercase):", transcriptLower);
      
      const mondayMatches = transcriptLower.match(/(hey monday|hey monty)/g);
      console.log("Monday matches found:", mondayMatches);
      
      if (mondayMatches && mondayMatches.length > 0) {
        const currentTime = Date.now();
        console.log("Monday trigger detected! Current time:", currentTime, "Last trigger time:", lastTriggerTime);
        
        if (currentTime - lastTriggerTime > 2000) { 
          console.log("Trigger cooldown passed, activating Monday...");
          setLastTriggerTime(currentTime);
          handleMondayActivation(transcript);
        } else {
          console.log("Trigger too recent, ignoring...");
        }
      } else {
        console.log("No Monday trigger found in transcript");
      }
    } else {
      // Debug why transcript processing is not happening
      if (!transcript) console.log("No transcript available");
      if (!listening) console.log("Not listening");
      if (conversationMode) console.log("In conversation mode, not looking for trigger");
      if (!isListeningForTrigger) console.log("Not listening for trigger");
    }
  }, [transcript, listening, conversationMode, isListeningForTrigger, lastTriggerTime, handleMondayActivation]);

  // Continuous Conversation Listening Effect
  useEffect(() => {
    if (conversationMode && !listening && isMounted && browserSupportsSpeechRecognition && isMicrophoneAvailable === true) {
      console.log("Starting listening for continuous conversation.");
      SpeechRecognition.startListening({ continuous: true, language: 'en-US', interimResults: true });
    }
  }, [conversationMode, listening, isMounted, browserSupportsSpeechRecognition, isMicrophoneAvailable]);

  // CRITICAL: Handle continuous conversation transcript processing
  useEffect(() => {
    if (conversationMode && transcript && listening) {
      setLiveTranscript(transcript);
      
      // Check for "Hey Monday" interruption ONLY - this resets to greeting
      const transcriptLower = transcript.toLowerCase().trim();
      if (transcriptLower.includes("hey monday") || transcriptLower.includes("hey monty")) {
        console.log("Hey Monday interruption detected - resetting conversation");
        if (isSpeaking) stopSpeaking();
        handleMondayActivation(transcript);
        return;
      }
      
      // Check for YouTube control commands
      const youtubeCommands = [
        "play video", "play", "pause video", "pause", "stop video", "stop",
        "mute video", "mute", "unmute video", "unmute", "restart video", "restart"
      ];
      
      if (youtubeVideoId && youtubeCommands.includes(transcriptLower)) {
        console.log("YouTube command detected:", transcriptLower);
        handleYouTubeCommand(transcriptLower);
        resetTranscript();
        setLiveTranscript("");
        return;
      }
      
      // Process complete phrases for conversation
      const trimmedTranscript = transcript.trim();
      if (trimmedTranscript.length > 0 && !isProcessingVoice) {
        const endsWithPunctuation = /[.!?]$/.test(trimmedTranscript);
        
        if (endsWithPunctuation || trimmedTranscript.includes('.') || trimmedTranscript.includes('?') || trimmedTranscript.includes('!')) {
          if (trimmedTranscript !== lastProcessedTranscript) {
            console.log("Processing sentence with punctuation:", trimmedTranscript);
            handleContinuousConversation(trimmedTranscript);
            resetTranscript();
          }
        } 
        else if (trimmedTranscript.length > 8 && trimmedTranscript !== lastProcessedTranscript) {
          const timeoutId = setTimeout(() => {
            if (trimmedTranscript === transcript.trim() && !isProcessingVoice) {
              console.log("Processing phrase after timeout:", trimmedTranscript);
              handleContinuousConversation(trimmedTranscript);
              resetTranscript();
            }
          }, 1500);
          
          return () => clearTimeout(timeoutId);
        }
      }
    }
  }, [conversationMode, transcript, listening, isProcessingVoice, handleContinuousConversation, lastProcessedTranscript, resetTranscript, youtubeVideoId, handleYouTubeCommand, isSpeaking, stopSpeaking, handleMondayActivation]);

  // WebSocket Message Handling Effect
  useEffect(() => {
    if (wsLastMessage?.data) {
      try {
        const message: WebSocketMessage = JSON.parse(wsLastMessage.data as string);
        switch (message.type) {
          case "voice_response": handleVoiceResponse(message); break;
          case "reasoning_progress": case "research_progress": handleProgressUpdate(message); break;
          case "voice_error": 
            setApiError(message.error || "An unknown voice error occurred");
            setIsThinking(false); setIsReceiving(false); setIsProcessingVoice(false);
            break;
          case "spatial_response": break;
          case "heartbeat_ack": break;
          default: console.warn("Unknown WebSocket message type received:", message.type);
        }
      } catch (e) {
        console.error("Failed to parse WebSocket message:", e, "Raw data:", wsLastMessage.data);
        setApiError("Invalid message received from server.");
      }
    }
  }, [wsLastMessage, handleVoiceResponse, handleProgressUpdate]); 

  // Typewriter Effect for Main Response
  useEffect(() => {
    if (isReceiving && fullResponseText && charIndex < fullResponseText.length) {
      const baseDelay = charIndex < 20 ? 2 : (charIndex < 50 ? 3 : 5);
      const randomDelay = charIndex < 50 ? 1 : 3;
      const timeout = setTimeout(() => {
        setCurrentResponse(prev => fullResponseText.slice(0, prev.length + 1));
        setCharIndex(prev => prev + 1);
      }, baseDelay + Math.random() * randomDelay);
      return () => clearTimeout(timeout);
    } else if (fullResponseText && charIndex >= fullResponseText.length && isReceiving) {
      setIsReceiving(false);
      setCurrentResponse(fullResponseText);
    }
  }, [charIndex, fullResponseText, isReceiving]);

  // --- Helper Functions ---
  const getModeConfig = (mode: string) => {
    const configs = {
      "reasoning": { icon: <Brain className="w-4 h-4" />, label: "Reasoning Pro", color: "from-purple-500 to-purple-700", description: "Deep analytical thinking" },
      "deep-research": { icon: <Activity className="w-4 h-4" />, label: "Deep Research", color: "from-blue-500 to-blue-700", description: "Comprehensive investigation" },
      "basic": { icon: <Search className="w-4 h-4" />, label: "Basic", color: "from-green-500 to-green-700", description: "Quick search & answers" }
    };
    return configs[mode as keyof typeof configs] || configs.basic;
  };
  
  // --- Render (JSX) ---
  return (
    <div className="min-h-screen bg-gradient-to-br from-[#091717] via-[#0A1A1A] to-[#091717] p-8 text-[#FBFAF4]">
      <div className="w-full max-w-[2000px] mx-auto">
        {/* Header and Connection Status */}
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
          <div className="flex items-center justify-center gap-2 mb-4">
            <div className={`w-3 h-3 rounded-full transition-colors duration-300 ${wsIsConnected ? 'bg-green-400 animate-pulse' : 'bg-red-400'}`}></div>
            <span className="text-sm text-[#20808D]">
              {wsIsConnected ? 'Connected to Monday Backend' : 
               apiError ? apiError : 
               wsHookError ? `Connection Error: ${ (wsHookError instanceof Error) ? wsHookError.message : (typeof wsHookError === 'string' ? wsHookError : 'General connection failure') }` : 
               'Connecting to Monday Backend...'}
            </span>
            {(!wsIsConnected && isMounted && (apiError || wsHookError)) && (
              <button 
                onClick={() => { 
                  console.log("FRONTEND: Manual reconnect initiated."); 
                  setApiError("Attempting to reconnect..."); 
                  if (wsReconnect) wsReconnect(); 
                  else console.warn("wsReconnect function not available from useMondayWebSocket hook.");
                }}
                className="ml-2 px-3 py-1 text-xs bg-[#20808D]/20 hover:bg-[#20808D]/30 border border-[#20808D]/40 rounded-lg transition-colors"
              >
                Try Reconnect
              </button>
            )}
            <div className={`w-3 h-3 rounded-full transition-colors duration-300 ${isInitialized ? 'bg-blue-400' : 'bg-gray-400'} ml-4`}></div>
            <span className="text-sm text-blue-400">{isInitialized ? 'Voice Ready' : 'Initializing Voice...'}</span>
          </div>
        </div>

        {/* Mode Buttons */}
        <div className="flex justify-center gap-4 mb-12">
          {(["basic", "reasoning", "deep-research"] as const).map((modeName) => {
            const config = getModeConfig(modeName);
            const isActive = currentMode === modeName;
            return (
              <button
                key={modeName}
                onClick={() => startNewResponseFlow(modeName)}
                className={`relative flex items-center gap-3 px-6 py-3 rounded-2xl transition-all duration-300 text-sm font-bold backdrop-blur-xl border transform hover:scale-105 ${isActive ? `bg-gradient-to-r ${config.color} text-white border-white/50 shadow-xl shadow-current/30 scale-110` : "bg-gradient-to-r from-[#20808D]/20 to-[#20808D]/10 text-[#20808D] hover:from-[#20808D]/30 hover:to-[#20808D]/20 border-[#20808D]/30"}`}
              >
                {isActive && isModeSwitching && <div className="absolute inset-0 rounded-2xl bg-gradient-to-r from-white/20 to-transparent animate-pulse"></div>}
                {config.icon}
                <div className="flex flex-col items-start">
                  <span>Sonar {config.label}</span>
                  {isActive && (<span className="text-xs opacity-80">{config.description}</span>)}
                </div>
                {isActive && (<div className="absolute -top-1 -right-1 w-3 h-3 bg-white rounded-full animate-pulse"></div>)}
              </button>
            );
          })}
        </div>

        {/* Mode Switch Message */}
        {modeSwitchMessage && (
          <div className="text-center mb-4 p-3 bg-gradient-to-r from-blue-500/20 to-purple-500/20 border border-blue-500/30 rounded-xl">
            <div className="flex items-center justify-center gap-2">
              <div className="w-4 h-4 border-2 border-blue-400 border-t-transparent rounded-full animate-spin"></div>
              <span className="text-blue-400 font-semibold">{modeSwitchMessage}</span>
            </div>
          </div>
        )}

        {/* Status Indicators */}
        {(isSpeaking || (isProcessingVoice && currentMode !== 'basic') || isModeSwitching) && (
          <div className="text-center mb-4 text-[#20808D] text-sm">
            {isSpeaking && <span className="mr-4">🔊 Monday is speaking...</span>}
            {isProcessingVoice && currentMode !== 'basic' && <span className="mr-4">{currentMode === 'reasoning' ? '🧠 Analyzing...' : '🔍 Researching...'}</span>}
            {isModeSwitching && <span>⚡ Switching modes...</span>}
          </div>
        )}

        {/* Voice Control Button and Status */}
        <div className="flex justify-center items-center gap-4 mb-8">
          <button 
            onClick={toggleListening}
            disabled={!isMounted || !browserSupportsSpeechRecognition || !isMicrophoneAvailable} 
            className={`p-4 rounded-full transition-colors duration-200 ease-in-out shadow-xl group focus:outline-none focus:ring-2 focus:ring-[#20808D] focus:ring-opacity-50 ${(listening || conversationMode) ? 'bg-green-500 hover:bg-green-600' : 'bg-gray-500 hover:bg-gray-600'} ${(!isMounted || !browserSupportsSpeechRecognition || !isMicrophoneAvailable) ? 'opacity-50 cursor-not-allowed' : ''}`}
            title={(listening || conversationMode) ? "Stop Listening" : (!isMicrophoneAvailable ? "Grant Mic Permission" : "Start Listening")}
          >
            {(listening || conversationMode) ? <Mic className="w-8 h-8 text-white animate-pulse" /> : <MicOff className="w-8 h-8 text-white" />}
          </button>
          
          {/* Audio initialization button */}
          {!isInitialized && (
            <button
              onClick={() => {
                console.log("Manual audio initialization triggered");
                initializeAudioContext();
                // Test audio with a simple beep
                try {
                  const context = new (window.AudioContext || (window as any).webkitAudioContext)();
                  const oscillator = context.createOscillator();
                  const gainNode = context.createGain();
                  oscillator.connect(gainNode);
                  gainNode.connect(context.destination);
                  gainNode.gain.value = 0.1;
                  oscillator.frequency.value = 440;
                  oscillator.start();
                  oscillator.stop(context.currentTime + 0.1);
                } catch (e) {
                  console.error("Test beep failed:", e);
                }
              }}
              className="p-4 rounded-full bg-blue-500 hover:bg-blue-600 transition-colors duration-200 shadow-xl"
              title="Initialize Audio"
            >
              <svg className="w-8 h-8 text-white" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15.536 8.464a5 5 0 010 7.072m2.828-9.9a9 9 0 010 12.728M5.586 15H4a1 1 0 01-1-1v-4a1 1 0 011-1h1.586l4.707-4.707C10.923 3.663 12 4.109 12 5v14c0 .891-1.077 1.337-1.707.707L5.586 15z" />
              </svg>
            </button>
          )}
          
          <div className="text-center">
            <div className="text-sm text-[#20808D] mb-1">
              {!isMounted ? "Voice system loading..." : 
               !browserSupportsSpeechRecognition ? "❌ Voice not supported in browser" : 
               !isMicrophoneAvailable ? "🎤 Click 'Grant Microphone Access' below" : 
               conversationMode ? `💬 In conversation (${getModeConfig(currentMode).label})` : 
               listening ? "🎤 Listening for 'Hey Monday' (hands-free)" : "🔄 Starting voice activation..."}
            </div>
            {isMounted && browserSupportsSpeechRecognition && !isMicrophoneAvailable && (
              <button
                onClick={async () => {
                  try {
                    console.log("Manual microphone permission request initiated...");
                    setApiError("Requesting microphone permission...");
                    const stream = await navigator.mediaDevices.getUserMedia({ audio: true });
                    console.log("Microphone permission granted successfully!");
                    stream.getTracks().forEach(track => track.stop()); 
                    setApiError(null);
                    // Also initialize audio context
                    initializeAudioContext();
                    console.log("Audio context initialized after permission grant");
                  } catch (err) {
                     console.error("Manual microphone permission denied:", err);
                     setApiError("Microphone permission denied. Please check browser settings and try again.");
                  }
                }}
                className="text-xs text-blue-400 underline hover:text-blue-300 mt-1"
              >
                Grant Microphone Access
              </button>
            )}
            
            {/* Voice Test Button */}
            {isMounted && browserSupportsSpeechRecognition && isMicrophoneAvailable && (
              <button
                onClick={() => {
                  console.log("=== VOICE TEST INITIATED ===");
                  console.log("Current transcript:", transcript);
                  console.log("Is listening:", listening);
                  console.log("Conversation mode:", conversationMode);
                  console.log("Processing voice:", isProcessingVoice);
                  
                  // Test speech recognition
                  if (!listening) {
                    console.log("Starting test listening...");
                    SpeechRecognition.startListening({ 
                      continuous: true, 
                      language: 'en-US', 
                      interimResults: true 
                    });
                  } else {
                    console.log("Stopping current listening...");
                    SpeechRecognition.stopListening();
                  }
                }}
                className="text-xs text-green-400 underline hover:text-green-300 mt-1 ml-2"
              >
                {listening ? "Stop Voice Test" : "Test Voice Recognition"}
              </button>
            )}
             
            {!conversationMode && !listening && isMicrophoneAvailable && (<div className="text-xs text-gray-400">Activating hands-free voice control...</div>)}
            {!conversationMode && listening && (<div className="text-xs text-gray-400">Say "Hey Monday" followed by your command (hands-free mode active)</div>)}
            {conversationMode && (<div className="text-xs text-gray-400">Continue speaking naturally or say "goodbye" to end</div>)}
          </div>
        </div>

        {/* Live Transcript Displays */}
        {isMounted && listening && transcript && !conversationMode && isListeningForTrigger && (
          <div className="text-center mb-4 p-4 bg-gradient-to-r from-green-500/10 to-green-600/10 border border-green-500/30 rounded-xl animate-pulse">
            <div className="flex items-center justify-center gap-2 mb-2">
              <Mic className="w-4 h-4 text-green-400 animate-pulse" />
              <span className="text-green-400 font-semibold text-sm">🎤 Voice Detected (Listening for Trigger)</span>
            </div>
            <p className="text-green-300/90 text-sm font-mono italic">"{transcript}"</p>
          </div>
        )}
        {isMounted && listening && liveTranscript && conversationMode && (
          <div className="text-center mb-4 p-4 bg-gradient-to-r from-green-500/20 to-green-600/20 border border-green-500/40 rounded-xl shadow-lg shadow-green-500/20">
            <div className="flex items-center justify-center gap-2 mb-2">
              <div className="relative"><Mic className="w-5 h-5 text-green-400 animate-pulse" /><div className="absolute -inset-1 bg-green-400 rounded-full opacity-20 animate-ping"></div></div>
              <span className="text-green-400 font-bold text-sm">LIVE RECOGNITION (Conversation)</span>
            </div>
            <p className="text-green-300 text-base font-mono font-semibold">"{liveTranscript}"</p>
          </div>
        )}

        {/* Error Display */}
        {apiError && !wsIsConnected && (
             <div className="text-center text-red-400 mb-4 p-4 bg-red-500/10 border border-red-500/30 rounded-xl">Error: {apiError}</div>
        )}
        
        {/* User Transcript Display */}
        {userTranscript && !isProcessingVoice && (
          <div className="text-center text-gray-400 mb-4 p-4 bg-[#20808D]/10 border border-[#20808D]/20 rounded-xl">
            <span className="font-semibold text-gray-300">Command processed:</span> {userTranscript}
            <div className="text-xs text-[#20808D]/70 mt-1">Mode: {getModeConfig(currentMode).label}</div>
          </div>
        )}

        {/* Main Content Panels */}
        <div className="grid grid-cols-1 xl:grid-cols-12 gap-8 h-[800px]">
          <div className="xl:col-span-4 transform xl:-rotate-2 xl:scale-95 origin-right">
            <div className="h-full shadow-2xl shadow-[#20808D]/20">
              <AdaptiveVisualizationPanel mode={currentMode} query={userTranscript} visualizationData={visualizationData} isThinking={isThinking} progressUpdates={progressUpdates} sources={progressSources} reasoning={progressReasoning} citations={progressCitations} thinking={visualizationData?.thinking} />
            </div>
          </div>
          <div className="xl:col-span-4 z-10">
            <div className="h-full shadow-2xl shadow-[#20808D]/30">
              <VoiceProcessingPanel currentResponse={currentResponse} isThinking={isThinking} isTyping={isReceiving} mode={currentMode} userTranscript={userTranscript} voiceActive={listening} processingError={apiError && wsIsConnected ? apiError : null} isModeSwitching={isModeSwitching} modeSwitchMessage={modeSwitchMessage} liveTranscript={liveTranscript} conversationMode={conversationMode} />
            </div>
          </div>
          <div className="xl:col-span-4 transform xl:rotate-2 xl:scale-95 origin-left">
            <div className="h-full shadow-2xl shadow-[#20808D]/20">
              <YouTubePanel 
                videoId={youtubeVideoId} 
                query={visualizationData?.query || userTranscript}
                mode={currentMode}
                responseContent={visualizationData?.content || currentResponse}
              />
            </div>
          </div>
        </div>

        {/* Footer */}
        <div className="text-center mt-12">
          <div className="inline-flex items-center gap-4 bg-gradient-to-r from-[#20808D]/10 to-[#20808D]/5 rounded-2xl px-8 py-4 backdrop-blur-xl border border-[#20808D]/20">
            <PerplexityLogo className="w-6 h-6" />
            <span className="text-[#20808D]/80 text-sm font-medium">Powered by Perplexity Sonar API • ElevenLabs TTS • WebXR Ready • Quest 2/3 Compatible</span>
          </div>
        </div>
      </div>
    </div>
  );
}

