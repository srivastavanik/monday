interface VoiceProfile {
  fundamentalFrequency: { min: number; max: number };
  formants: number[];
  spectralCentroid: number;
  voiceId: string;
}

interface AnalysisResult {
  isMondayVoice: boolean;
  confidence: number;
  features: any;
}

class VoiceSignatureFilter {
  private mondayVoiceProfile: VoiceProfile;
  private mondayResponsePatterns: RegExp[] = [];
  
  constructor() {
    // Monday's voice characteristics from ElevenLabs
    this.mondayVoiceProfile = {
      fundamentalFrequency: { min: 80, max: 180 }, // Male voice range
      formants: [700, 1100, 2900], // Typical formant frequencies
      spectralCentroid: 1500, // Brightness characteristic
      voiceId: 'monday_elevenlabs_voice'
    };
    
    this.initializeResponsePatterns();
  }
  
  private initializeResponsePatterns(): void {
    // Patterns that indicate AI-generated responses
    this.mondayResponsePatterns = [
      /^(hello|hi|greetings),?\s+(i'm|i am)\s+monday/i,
      /^(that's|this is)\s+(a\s+)?(great|excellent|interesting)\s+question/i,
      /^(let me|i'll|i can|i will)\s+(explain|tell you|share|help)/i,
      /^(in\s+)?(quantum|physics|science|the field of)/i,
      /^(according to|based on|research shows|studies indicate)/i,
      /^(would you like|do you want|are you interested in)\s+(to know|learning)/i,
      /^(quantum mechanics|quantum physics|quantum computing)/i,
      /^(the answer is|the solution involves|this relates to)/i,
      /^(fascinating|interesting|remarkable)\s+(question|topic|concept)/i,
      /^(i found|i discovered|i learned)\s+(some|great|interesting)/i
    ];
  }
  
  // Quick text-based filtering for obvious AI responses
  isLikelyMondayResponse(transcript: string): boolean {
    const cleaned = transcript.trim().toLowerCase();
    
    // Check against known Monday response patterns
    for (const pattern of this.mondayResponsePatterns) {
      if (pattern.test(cleaned)) {
        console.log('VoiceSignatureFilter: 🚫 Detected Monday response pattern:', pattern.source);
        return true;
      }
    }
    
    // Check for characteristic Monday phrases - ONLY ACTUAL AI RESPONSE PHRASES
    const mondayPhrases = [
      "i'm monday",
      "hello, i'm monday", 
      "monday here",
      "let me share what i discovered",
      "i found some great information",
      "based on my research",
      "according to my analysis",
      "i've gathered more details",
      "that's a great topic",
      "thanks for sharing your curiosity"
    ];
    
    for (const phrase of mondayPhrases) {
      if (cleaned.includes(phrase)) {
        console.log('VoiceSignatureFilter: 🚫 Detected Monday phrase:', phrase);
        return true;
      }
    }
    
    // Check for AI-like response structure
    if (this.hasAIResponseStructure(cleaned)) {
      console.log('VoiceSignatureFilter: 🚫 Detected AI response structure');
      return true;
    }
    
    return false;
  }
  
  private hasAIResponseStructure(text: string): boolean {
    // AI responses often have these characteristics:
    
    // 1. Very formal/academic language - but only if combined with other indicators
    const formalWords = ['furthermore', 'moreover', 'consequently', 'therefore', 'specifically', 'particularly'];
    const formalCount = formalWords.filter(word => text.includes(word)).length;
    
    // 2. Multiple technical terms in sequence - RELAXED: Only flag if it's clearly an explanation
    const technicalTerms = ['quantum', 'mechanics', 'physics', 'theory', 'principle', 'phenomenon', 'research', 'study', 'analysis'];
    const technicalCount = technicalTerms.filter(term => text.includes(term)).length;
    
    // 3. Very long, complex sentences (typical of AI explanations)
    const sentences = text.split(/[.!?]+/).filter(s => s.trim().length > 0);
    const avgLength = sentences.reduce((sum, s) => sum + s.length, 0) / sentences.length;
    const isVeryLong = avgLength > 150 && sentences.length > 3; // Increased threshold
    
    // 4. Starts with explanatory phrases
    const explanatoryStarts = [
      'in essence', 'essentially', 'fundamentally', 'basically', 'in simple terms',
      'to put it simply', 'in other words', 'what this means is', 'this refers to'
    ];
    const hasExplanatoryStart = explanatoryStarts.some(start => text.startsWith(start));
    
    // 5. Contains AI-like response indicators
    const aiIndicators = [
      'let me explain', 'i can help you understand', 'here\'s what you need to know',
      'according to my analysis', 'based on my research', 'i\'ve found that'
    ];
    const hasAIIndicator = aiIndicators.some(indicator => text.includes(indicator));
    
    // CRITICAL FIX: Only flag as AI response if multiple indicators are present
    // Don't flag simple user commands that happen to contain technical terms
    
    // Short commands (< 50 chars) are likely user commands, not AI responses
    if (text.length < 50) return false;
    
    // User commands typically start with action words
    const userCommandStarts = ['can you', 'please', 'tell me', 'explain', 'what is', 'how does', 'research', 'think about'];
    const isUserCommand = userCommandStarts.some(start => text.startsWith(start));
    if (isUserCommand) return false;
    
    // Only flag if we have multiple strong indicators
    let indicatorCount = 0;
    if (formalCount >= 3) indicatorCount++; // Increased threshold
    if (technicalCount >= 5 && isVeryLong) indicatorCount++; // Only if very technical AND very long
    if (hasExplanatoryStart) indicatorCount++;
    if (hasAIIndicator) indicatorCount++;
    
    // Need at least 2 strong indicators to flag as AI response
    return indicatorCount >= 2;
  }
  
  // Enhanced recognition with filtering
  createFilteredRecognition(): any {
    // This would be implemented to wrap the native recognition
    // For now, we'll integrate this into the existing VoiceSystemController
    return null;
  }
  
  // Analyze audio characteristics (placeholder for future implementation)
  async analyzeAudioChunk(audioData: Float32Array): Promise<AnalysisResult> {
    // This would implement actual audio analysis
    // For now, return a basic result
    return {
      isMondayVoice: false,
      confidence: 0.5,
      features: {}
    };
  }
  
  // Get filter statistics
  getFilterStats(): any {
    return {
      patternsCount: this.mondayResponsePatterns.length,
      voiceProfile: this.mondayVoiceProfile
    };
  }
}

export { VoiceSignatureFilter };
export type { VoiceProfile, AnalysisResult }; 