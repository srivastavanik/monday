interface ConversationEntry {
  role: 'user' | 'assistant';
  content: string;
  ttsContent?: string;
  timestamp: number;
}

export class ContextCleaner {
  // Remove voice recognition artifacts and clean messages
  static cleanMessage(message: string): string {
    // Remove common artifacts
    const artifacts = [
      /^s for you to explore$/i,
      /^okay$/i,
      /^mm-hmm$/i,
      /^uh$/i,
      /^\w{1,2}$/i, // Single or double character fragments
      /^thanks for sharing your curiosity/i, // Monday's own phrases being picked up
      /^i've gathered more details/i,
      /^that's a great topic/i,
    ];
    
    for (const artifact of artifacts) {
      if (artifact.test(message.trim())) {
        return ''; // Return empty to be filtered out
      }
    }
    
    // Clean up truncated sentences
    if (message.endsWith('...') || message.endsWith('…')) {
      // This is fine, keep it
    } else if (message.length < 10 && !message.endsWith('.') && !message.endsWith('?') && !message.endsWith('!')) {
      // Likely a fragment
      return '';
    }
    
    return message.trim();
  }
  
  static validateAndCleanContext(entries: ConversationEntry[]): Array<{role: 'user' | 'assistant', content: string}> {
    const cleaned: Array<{role: 'user' | 'assistant', content: string}> = [];
    
    for (const entry of entries) {
      const cleanContent = this.cleanMessage(entry.content);
      
      // Skip empty or invalid messages
      if (!cleanContent || cleanContent.length < 3) {
        continue;
      }
      
      cleaned.push({
        role: entry.role,
        content: cleanContent
      });
    }
    
    return cleaned;
  }
  
  // Convert old string format to new format for backward compatibility
  static convertLegacyContext(legacyContext: string[]): ConversationEntry[] {
    const entries: ConversationEntry[] = [];
    
    for (const ctx of legacyContext) {
      if (ctx.startsWith('User: ')) {
        entries.push({
          role: 'user',
          content: ctx.replace('User: ', ''),
          timestamp: Date.now()
        });
      } else if (ctx.startsWith('Monday: ')) {
        entries.push({
          role: 'assistant',
          content: ctx.replace('Monday: ', ''),
          timestamp: Date.now()
        });
      }
    }
    
    return entries;
  }
} 