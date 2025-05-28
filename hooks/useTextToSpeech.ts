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
          console.log('Audio context resume failed:', error);
        });
      }
    } catch (error) {
      console.log('Audio context initialization failed:', error);
    }
  }, []);

  const speak = useCallback(async (text: string) => {
    console.log("TTS speak called with:", { 
      textLength: text?.length, 
      isInitialized, 
      isSpeaking, 
      hasApiKey: !!apiKey,
      voiceId 
    });
    
    if (!isInitialized) {
      console.log('Audio context not initialized, skipping TTS');
      return;
    }

    if (isSpeaking) {
      console.log('Already speaking, ignoring new request');
      return;
    }

    if (!apiKey || !voiceId) {
      console.error('Missing ElevenLabs API key or voice ID, cannot speak');
      return;
    }

    setIsSpeaking(true);
    
    try {
      console.log('Making request to ElevenLabs API...');
      
      const response = await fetch(`https://api.elevenlabs.io/v1/text-to-speech/${voiceId}`, {
        method: 'POST',
        headers: {
          'Accept': 'audio/mpeg',
          'Content-Type': 'application/json',
          'xi-api-key': apiKey
        },
        body: JSON.stringify({
          text,
          model_id: 'eleven_monolingual_v1',
          voice_settings: {
            stability: 0.5,
            similarity_boost: 0.5
          }
        })
      });

      if (!response.ok) {
        if (response.status === 401) {
          console.error('ElevenLabs API authentication failed (401). Please check your API key.');
          console.error('The API key may be invalid, expired, or incorrectly formatted.');
          console.error('To get a valid API key:');
          console.error('1. Sign up or log in at https://elevenlabs.io');
          console.error('2. Go to your profile settings');
          console.error('3. Navigate to the API Keys section');
          console.error('4. Generate a new API key and update it in the code');
          console.warn('Continuing without voice output...');
        } else {
          console.error(`ElevenLabs API error: ${response.status} ${response.statusText}`);
        }
        throw new Error(`ElevenLabs API error: ${response.status}`);
      }

      console.log('ElevenLabs API response received, creating audio...');
      
      const audioBlob = await response.blob();
      const audioUrl = URL.createObjectURL(audioBlob);
      
      const audio = new Audio(audioUrl);
      currentAudioRef.current = audio;
      
      audio.volume = 0.8;
      
      audio.onended = () => {
        console.log('ElevenLabs audio playback ended');
        setIsSpeaking(false);
        URL.revokeObjectURL(audioUrl);
        currentAudioRef.current = null;
      };
      
      audio.onerror = (e) => {
        console.error('ElevenLabs audio playback error:', e);
        setIsSpeaking(false);
        URL.revokeObjectURL(audioUrl);
        currentAudioRef.current = null;
      };
      
      console.log('Starting ElevenLabs audio playback...');
      await audio.play();
      console.log('ElevenLabs audio playback started successfully');
      
    } catch (error) {
      console.error('ElevenLabs TTS error:', error);
      setIsSpeaking(false);
      // Don't throw the error, just log it so the app continues to work
    }
  }, [voiceId, apiKey, isSpeaking, isInitialized]);

  const speakProgress = useCallback(async (text: string) => {
    // For progress updates, we'll use a quieter, shorter version
    if (!isInitialized || isSpeaking || !apiKey || !voiceId) {
      return;
    }

    try {
      // Use ElevenLabs for progress updates too
      const response = await fetch(`https://api.elevenlabs.io/v1/text-to-speech/${voiceId}`, {
        method: 'POST',
        headers: {
          'Accept': 'audio/mpeg',
          'Content-Type': 'application/json',
          'xi-api-key': apiKey
        },
        body: JSON.stringify({
          text,
          model_id: 'eleven_monolingual_v1',
          voice_settings: {
            stability: 0.7,
            similarity_boost: 0.3
          }
        })
      });

      if (!response.ok) {
        if (response.status === 401) {
          console.error('ElevenLabs API authentication failed for progress update');
          // Silently fail for progress updates to avoid spamming console
        }
        throw new Error(`ElevenLabs API error: ${response.status}`);
      }

      const audioBlob = await response.blob();
      const audioUrl = URL.createObjectURL(audioBlob);
      
      // Create new audio element for progress
      const audio = new Audio(audioUrl);
      audio.volume = 0.5; // Quieter for progress updates
      
      audio.onended = () => {
        URL.revokeObjectURL(audioUrl);
      };
      
      await audio.play();
    } catch (error) {
      console.error('Progress TTS error:', error);
      // No fallback - just log the error
    }
  }, [voiceId, apiKey, isSpeaking, isInitialized]);

  const stop = useCallback(() => {
    try {
      // Stop ElevenLabs audio
      if (currentAudioRef.current) {
        currentAudioRef.current.pause();
        currentAudioRef.current = null;
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