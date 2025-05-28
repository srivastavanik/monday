"use client"

import { useState, useCallback, useRef } from 'react';

interface UseTextToSpeechProps {
  voiceId?: string;
  apiKey?: string;
}

export interface UseTextToSpeechReturn {
  speak: (text: string) => Promise<void>;
  speakProgress: (text: string) => Promise<void>;
  stop: () => void;
  isSpeaking: boolean;
  initializeAudioContext: () => void;
  isInitialized: boolean;
}

export function useTextToSpeech({ voiceId, apiKey }: UseTextToSpeechProps = {}) {
  const [isSpeaking, setIsSpeaking] = useState(false);
  const [isInitialized, setIsInitialized] = useState(true); // Always initialized for hands-free
  const audioContextRef = useRef<AudioContext | null>(null);
  const currentAudioRef = useRef<HTMLAudioElement | null>(null);

  const initializeAudioContext = useCallback(() => {
    if (typeof window === 'undefined') return;
    
    try {
      // Always set as initialized for hands-free operation
      setIsInitialized(true);
      
      // Try to initialize audio context but don't block if it fails
      if (!audioContextRef.current) {
        audioContextRef.current = new (window.AudioContext || (window as any).webkitAudioContext)();
      }
      
      // Try to resume audio context but don't block on it
      if (audioContextRef.current.state === 'suspended') {
        audioContextRef.current.resume().catch((error) => {
          console.log('Audio context resume failed, using Web Speech API:', error);
        });
      }
    } catch (error) {
      console.log('Audio context initialization failed, using Web Speech API:', error);
    }
  }, []);

  const fallbackToWebSpeech = useCallback((text: string) => {
    if (typeof window === 'undefined' || !('speechSynthesis' in window)) {
      console.log('Web Speech API not available');
      setIsSpeaking(false);
      return;
    }

    try {
      // Cancel any ongoing speech
      window.speechSynthesis.cancel();
      
      const utterance = new SpeechSynthesisUtterance(text);
      utterance.rate = 0.9;
      utterance.pitch = 1;
      utterance.volume = 0.8;
      
      utterance.onstart = () => {
        console.log('Web Speech API started speaking');
      };
      
      utterance.onend = () => {
        console.log('Web Speech API finished speaking');
        setIsSpeaking(false);
      };
      
      utterance.onerror = (event) => {
        console.log('Web Speech API error:', event.error);
        setIsSpeaking(false);
      };
      
      // Speak immediately for hands-free operation
      window.speechSynthesis.speak(utterance);
      
    } catch (error) {
      console.log('Web Speech API failed:', error);
      setIsSpeaking(false);
    }
  }, []);

  const speak = useCallback(async (text: string) => {
    if (!text.trim()) return;
    
    try {
      // Stop any current speech
      if (currentAudioRef.current) {
        currentAudioRef.current.pause();
        currentAudioRef.current = null;
      }
      
      setIsSpeaking(true);
      
      // Always try ElevenLabs first if credentials are provided
      if (apiKey && voiceId) {
        try {
          console.log('Attempting to use ElevenLabs TTS...');
          const response = await fetch(`https://api.elevenlabs.io/v1/text-to-speech/${voiceId}`, {
            method: 'POST',
            headers: {
              'Accept': 'audio/mpeg',
              'Content-Type': 'application/json',
              'xi-api-key': apiKey,
            },
            body: JSON.stringify({
              text: text,
              model_id: 'eleven_monolingual_v1',
              voice_settings: {
                stability: 0.5,
                similarity_boost: 0.5,
              },
            }),
          });

          if (response.ok) {
            console.log('ElevenLabs TTS response received successfully');
            const audioBlob = await response.blob();
            const audioUrl = URL.createObjectURL(audioBlob);
            const audio = new Audio(audioUrl);
            currentAudioRef.current = audio;
            
            // Set up promise-based playback
            return new Promise<void>((resolve, reject) => {
              audio.onended = () => {
                console.log('ElevenLabs audio finished playing');
                setIsSpeaking(false);
                URL.revokeObjectURL(audioUrl);
                currentAudioRef.current = null;
                resolve();
              };
              
              audio.onerror = (error) => {
                console.log('ElevenLabs audio playback error:', error);
                setIsSpeaking(false);
                currentAudioRef.current = null;
                URL.revokeObjectURL(audioUrl);
                reject(error);
              };
              
              // Try to play the audio
              console.log('Playing ElevenLabs audio...');
              const playPromise = audio.play();
              
              if (playPromise !== undefined) {
                playPromise
                  .then(() => {
                    console.log('ElevenLabs audio started playing successfully');
                  })
                  .catch((error) => {
                    console.log('ElevenLabs play() promise rejected:', error);
                    // Don't immediately fall back - let the audio element try to recover
                    // Only reject if it's a real error, not a temporary interruption
                    if (error.name !== 'AbortError') {
                      reject(error);
                    }
                  });
              }
            });
          } else {
            console.log(`ElevenLabs API error: ${response.status} ${response.statusText}`);
            const errorText = await response.text();
            console.log('ElevenLabs error response:', errorText);
            throw new Error(`ElevenLabs API error: ${response.status}`);
          }
        } catch (error) {
          console.log('ElevenLabs TTS failed, falling back to Web Speech API:', error);
          // Fall back to Web Speech API
          fallbackToWebSpeech(text);
        }
      } else {
        console.log('No ElevenLabs credentials provided, using Web Speech API');
        // Use Web Speech API if no ElevenLabs credentials
        fallbackToWebSpeech(text);
      }
      
    } catch (error) {
      console.error('Text-to-speech error:', error);
      setIsSpeaking(false);
    }
  }, [apiKey, voiceId, fallbackToWebSpeech]);

  // Speak progress updates (shorter, less intrusive)
  const speakProgress = useCallback(async (text: string) => {
    if (!isInitialized) {
      console.warn('Audio context not initialized. Skipping progress speech.');
      return;
    }

    // Skip if already speaking main content
    if (isSpeaking) {
      console.log('Already speaking, skipping progress update');
      return;
    }

    // Use a simpler, faster approach for progress updates
    try {
      // Try Web Speech API first for quick progress updates
      if ('speechSynthesis' in window) {
        const utterance = new SpeechSynthesisUtterance(text);
        utterance.rate = 1.2; // Slightly faster for progress
        utterance.pitch = 0.9; // Slightly lower pitch
        utterance.volume = 0.8; // Slightly quieter
        
        speechSynthesis.speak(utterance);
        return;
      }
    } catch (error) {
      console.log('Progress speech failed, continuing silently');
    }
  }, [isInitialized, isSpeaking]);

  const stop = useCallback(() => {
    try {
      // Stop ElevenLabs audio
      if (currentAudioRef.current) {
        currentAudioRef.current.pause();
        currentAudioRef.current = null;
      }
      
      // Stop Web Speech API
      if (typeof window !== 'undefined' && 'speechSynthesis' in window) {
        window.speechSynthesis.cancel();
      }
      
      setIsSpeaking(false);
    } catch (error) {
      console.log('Error stopping speech:', error);
      setIsSpeaking(false);
    }
  }, []);

  return {
    speak,
    speakProgress,
    stop,
    isSpeaking,
    initializeAudioContext,
    isInitialized
  };
} 