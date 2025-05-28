"use client"

import { useState, useCallback, useRef, useEffect } from 'react';

interface TextToSpeechOptions {
  voiceId?: string;
  modelId?: string;
  apiKey?: string;
}

export const useTextToSpeech = (options: TextToSpeechOptions = {}) => {
  const [isSpeaking, setIsSpeaking] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const audioRef = useRef<HTMLAudioElement | null>(null);
  const audioContextRef = useRef<AudioContext | null>(null);

  // Initialize AudioContext on user interaction
  const initializeAudioContext = useCallback(() => {
    if (!audioContextRef.current) {
      audioContextRef.current = new (window.AudioContext || (window as any).webkitAudioContext)();
    }
  }, []);

  // Use ElevenLabs API for text-to-speech
  const speak = useCallback(async (text: string) => {
    if (!text) return;

    setError(null);
    setIsLoading(true);

    try {
      const apiKey = process.env.NEXT_PUBLIC_ELEVENLABS_API_KEY || options.apiKey;
      const voiceId = options.voiceId || "21m00Tcm4TlvDq8ikWAM"; // Default to Rachel voice
      
      if (!apiKey) {
        throw new Error('ElevenLabs API key not found');
      }

      // Call ElevenLabs API
      const response = await fetch(`https://api.elevenlabs.io/v1/text-to-speech/${voiceId}`, {
        method: 'POST',
        headers: {
          'Accept': 'audio/mpeg',
          'Content-Type': 'application/json',
          'xi-api-key': apiKey,
        },
        body: JSON.stringify({
          text: text,
          model_id: options.modelId || "eleven_monolingual_v1",
          voice_settings: {
            stability: 0.5,
            similarity_boost: 0.75,
            style: 0.0,
            use_speaker_boost: true
          }
        }),
      });

      if (!response.ok) {
        throw new Error(`ElevenLabs API error: ${response.status}`);
      }

      const audioBlob = await response.blob();
      const audioUrl = URL.createObjectURL(audioBlob);
      
      // Create and play audio
      if (audioRef.current) {
        audioRef.current.pause();
        URL.revokeObjectURL(audioRef.current.src);
      }

      audioRef.current = new Audio(audioUrl);
      
      audioRef.current.onloadstart = () => {
        setIsSpeaking(true);
        setIsLoading(false);
      };

      audioRef.current.onended = () => {
        setIsSpeaking(false);
        URL.revokeObjectURL(audioUrl);
      };

      audioRef.current.onerror = (event) => {
        setError('Audio playback error');
        setIsSpeaking(false);
        setIsLoading(false);
        URL.revokeObjectURL(audioUrl);
      };

      await audioRef.current.play();
      
    } catch (err) {
      console.error('TTS Error:', err);
      setError(err instanceof Error ? err.message : 'Failed to synthesize speech');
      setIsLoading(false);
      setIsSpeaking(false);
    }
  }, [options.apiKey, options.voiceId, options.modelId]);

  const stop = useCallback(() => {
    if (audioRef.current) {
      audioRef.current.pause();
      audioRef.current.currentTime = 0;
      if (audioRef.current.src) {
        URL.revokeObjectURL(audioRef.current.src);
      }
    }
    setIsSpeaking(false);
  }, []);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      stop();
    };
  }, [stop]);

  return {
    speak,
    stop,
    isSpeaking,
    isLoading,
    error,
    initializeAudioContext
  };
}; 