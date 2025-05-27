interface TTSRecord {
  text: string;
  timestamp: number;
  audioFeatures?: AudioFeatures;
}

interface CommandRecord {
  text: string;
  timestamp: number;
  processed: boolean;
}

interface AudioFeatures {
  pitch?: number;
  spectralCentroid?: number;
  mfcc?: number[];
}

interface ValidationResult {
  valid: boolean;
  reason?: string;
  confidence: number;
}

interface Validation {
  valid: boolean;
  reason?: string;
  confidence?: number;
}

class CommandValidator {
  private recentTTSOutputs: TTSRecord[] = [];
  private commandHistory: CommandRecord[] = [];
  private feedbackPatterns: RegExp[] = [];
  private mondayPhrases: string[] = [];
  
  constructor() {
    this.initializeFeedbackPatterns();
    this.initializeMondayPhrases();
  }
  
  private initializeFeedbackPatterns(): void {
    // Common feedback loop patterns that indicate the system is hearing itself
    // Made more specific to avoid blocking legitimate user commands
    this.feedbackPatterns = [
      /^(okay|alright|sure|yes|mm-hmm|uh-huh)\s*$/i,
      /^(i understand|i see|got it|i get it)\s*$/i,
      /^(that's|this is|it's)\s+(correct|right|true|interesting)/i,
      /^(to|the|a|an)\s+\w+\s+(is|are|was|were)/i,
      /^(let me|i'll|i can|i will)\s+(explain|tell|show)/i,
      /^(would you like|do you want|are you interested)/i,
      // Removed overly broad patterns that block legitimate commands:
      // /^(quantum|mechanics|physics|science|research)/i, // TOO BROAD - blocks user requests
      /^(monday|hello|hi there|greetings)/i,
      /^(that's a great question|excellent question)/i,
      /^(in quantum mechanics|in physics|in science)/i,
      /^(the answer is|the solution is|the result is)/i,
      // More specific patterns for actual AI responses:
      /^(according to|based on)\s+(my research|the research|studies)/i,
      /^(i found some|i discovered some|here's what i)/i,
      /^(let me share what|i'll explain|allow me to)/i
    ];
  }
  
  private initializeMondayPhrases(): void {
    // Phrases that Monday commonly uses - strong indicators of feedback
    // Made more specific to avoid blocking legitimate user commands
    this.mondayPhrases = [
      "i'm monday",
      "hello, i'm monday", 
      "monday here",
      "let me explain",
      "that's a great question",
      "would you like to know more",
      "in the quantum world",
      "let me share what i discovered",
      "i found some great information",
      "here's what i learned",
      "based on my research",
      "according to the latest studies",
      // More specific Monday response patterns:
      "i've gathered more details",
      "hey there. i've gathered",
      "that's a great topic to think through",
      "while i'm having some connectivity issues",
      // Only keep very specific AI response phrases:
      "i can help you understand",
      "here's what you need to know",
      "according to my analysis",
      "i've found that",
      "thanks for sharing your curiosity"
    ];
  }
  
  async validateCommand(command: string, audioFeatures?: AudioFeatures): Promise<ValidationResult> {
    console.log('CommandValidator: 🔍 Validating command:', command);
    
    const validations = await Promise.all([
      this.checkAgainstRecentTTS(command),
      this.checkFeedbackPattern(command),
      this.checkMondayPhrases(command),
      this.checkTemporalPattern(command),
      this.checkCommandLength(command),
      this.checkRepetitivePattern(command)
    ]);
    
    const issues = validations.filter(v => !v.valid);
    
    if (issues.length > 0) {
      const primaryIssue = issues[0];
      console.log('CommandValidator: 🚫 Command validation failed:', primaryIssue.reason);
      return {
        valid: false,
        reason: primaryIssue.reason,
        confidence: primaryIssue.confidence || 0.8
      };
    }
    
    // Record valid command
    this.recordCommand(command);
    console.log('CommandValidator: ✅ Command validated successfully');
    return { valid: true, confidence: 1.0 };
  }
  
  private async checkAgainstRecentTTS(command: string): Promise<Validation> {
    // Check against recent TTS outputs with fuzzy matching
    for (const tts of this.recentTTSOutputs) {
      const similarity = this.calculateTextSimilarity(command, tts.text);
      if (similarity > 0.8) { // Increased threshold to 80% to be less aggressive
        return {
          valid: false,
          reason: `Command too similar to recent TTS output (${Math.round(similarity * 100)}% match)`,
          confidence: similarity
        };
      }
      
      // Check for partial matches - but only if it's a very close match and short
      if (tts.text.toLowerCase().includes(command.toLowerCase()) && 
          command.length > 15 && // Increased minimum length
          similarity > 0.7) { // Added similarity check
        return {
          valid: false,
          reason: 'Command appears to be fragment of recent TTS output',
          confidence: 0.9
        };
      }
    }
    return { valid: true };
  }
  
  private checkFeedbackPattern(command: string): Validation {
    // Check against known feedback patterns
    for (const pattern of this.feedbackPatterns) {
      if (pattern.test(command.trim())) {
        return {
          valid: false,
          reason: 'Command matches known feedback pattern',
          confidence: 0.9
        };
      }
    }
    return { valid: true };
  }
  
  private checkMondayPhrases(command: string): Validation {
    // Check for Monday's characteristic phrases
    const lowerCommand = command.toLowerCase();
    for (const phrase of this.mondayPhrases) {
      if (lowerCommand.includes(phrase)) {
        return {
          valid: false,
          reason: `Command contains Monday's characteristic phrase: "${phrase}"`,
          confidence: 0.95
        };
      }
    }
    return { valid: true };
  }
  
  private checkTemporalPattern(command: string): Validation {
    // Check for suspicious timing patterns
    const now = Date.now();
    const recentCommands = this.commandHistory.filter(
      c => now - c.timestamp < 10000 // Last 10 seconds
    );
    
    // Too many rapid commands (likely TTS chunks) - increased threshold
    if (recentCommands.length > 8) { // Increased from 5 to 8
      return {
        valid: false,
        reason: 'Too many rapid commands detected (likely feedback)',
        confidence: 0.8
      };
    }
    
    // Check for very recent command (within 1 second) - reduced from 2 seconds
    const veryRecentCommands = recentCommands.filter(
      c => now - c.timestamp < 1000 // Reduced to 1 second
    );
    
    if (veryRecentCommands.length > 2) { // Increased threshold
      return {
        valid: false,
        reason: 'Commands arriving too rapidly (likely feedback)',
        confidence: 0.85
      };
    }
    
    return { valid: true };
  }
  
  private checkCommandLength(command: string): Validation {
    // Very short commands during active conversation are suspicious
    if (command.trim().length < 3) {
      return {
        valid: false,
        reason: 'Command too short (likely audio artifact)',
        confidence: 0.7
      };
    }
    
    // Very long commands that look like TTS chunks
    if (command.length > 200) {
      const sentences = command.split(/[.!?]+/).length;
      if (sentences > 3) {
        return {
          valid: false,
          reason: 'Command too long and complex (likely TTS output)',
          confidence: 0.8
        };
      }
    }
    
    return { valid: true };
  }
  
  private checkRepetitivePattern(command: string): Validation {
    // Check if this exact command was recently processed
    const recentExact = this.commandHistory.find(
      c => c.text.toLowerCase() === command.toLowerCase() && 
           Date.now() - c.timestamp < 10000 // Reduced to 10 seconds - allow repeated questions after 10s
    );
    
    if (recentExact) {
      // Allow repeated legitimate user commands (thinking, research requests)
      const legitimateRepeats = [
        /^(can you )?think about/i,
        /^(can you )?research/i,
        /^(hey )?monday/i,
        /^(what|how|why|when|where)/i
      ];
      
      for (const pattern of legitimateRepeats) {
        if (pattern.test(command)) {
          console.log('CommandValidator: ✅ Allowing repeated legitimate command:', command);
          return { valid: true };
        }
      }
      
      return {
        valid: false,
        reason: 'Exact command recently processed (likely feedback loop)',
        confidence: 0.9
      };
    }
    
    return { valid: true };
  }
  
  private calculateTextSimilarity(text1: string, text2: string): number {
    // Simple Levenshtein distance-based similarity
    const longer = text1.length > text2.length ? text1 : text2;
    const shorter = text1.length > text2.length ? text2 : text1;
    
    if (longer.length === 0) return 1.0;
    
    const distance = this.levenshteinDistance(longer.toLowerCase(), shorter.toLowerCase());
    return (longer.length - distance) / longer.length;
  }
  
  private levenshteinDistance(str1: string, str2: string): number {
    const matrix = Array(str2.length + 1).fill(null).map(() => Array(str1.length + 1).fill(null));
    
    for (let i = 0; i <= str1.length; i++) matrix[0][i] = i;
    for (let j = 0; j <= str2.length; j++) matrix[j][0] = j;
    
    for (let j = 1; j <= str2.length; j++) {
      for (let i = 1; i <= str1.length; i++) {
        const indicator = str1[i - 1] === str2[j - 1] ? 0 : 1;
        matrix[j][i] = Math.min(
          matrix[j][i - 1] + 1,     // deletion
          matrix[j - 1][i] + 1,     // insertion
          matrix[j - 1][i - 1] + indicator // substitution
        );
      }
    }
    
    return matrix[str2.length][str1.length];
  }
  
  // Public methods for recording TTS and commands
  recordTTSOutput(text: string, audioFeatures?: AudioFeatures): void {
    this.recentTTSOutputs.push({
      text,
      timestamp: Date.now(),
      audioFeatures
    });
    
    // Keep only recent TTS outputs (last 2 minutes)
    const cutoff = Date.now() - 120000;
    this.recentTTSOutputs = this.recentTTSOutputs.filter(tts => tts.timestamp > cutoff);
    
    console.log('CommandValidator: 📝 Recorded TTS output for validation');
  }
  
  private recordCommand(command: string): void {
    this.commandHistory.push({
      text: command,
      timestamp: Date.now(),
      processed: true
    });
    
    // Keep only recent commands (last 5 minutes)
    const cutoff = Date.now() - 300000;
    this.commandHistory = this.commandHistory.filter(cmd => cmd.timestamp > cutoff);
  }
  
  // Get statistics for debugging
  getValidationStats(): any {
    return {
      recentTTSCount: this.recentTTSOutputs.length,
      recentCommandCount: this.commandHistory.length,
      lastTTS: this.recentTTSOutputs[this.recentTTSOutputs.length - 1]?.text.substring(0, 50),
      lastCommand: this.commandHistory[this.commandHistory.length - 1]?.text
    };
  }
  
  // Clear history (for testing or reset)
  clearHistory(): void {
    this.recentTTSOutputs = [];
    this.commandHistory = [];
    console.log('CommandValidator: 🧹 Cleared validation history');
  }
}

export { CommandValidator };
export type { ValidationResult, TTSRecord, CommandRecord, AudioFeatures }; 