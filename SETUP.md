# Monday Voice Agent - Final Setup

This is the working implementation of the Monday voice agent with all fixes applied.

## Quick Start

### 1. Navigate to Working Directory
```bash
cd nidsmonday
```

### 2. Install Dependencies
```bash
npm install
```

### 3. Start Backend (Terminal 1)
```bash
node backend-server.js
```
You should see: "Monday backend server running on port 3001"

### 4. Start Frontend (Terminal 2)
```bash
npm run dev
```
Go to: http://localhost:3000

### 5. Use Voice Commands
- Say: "Hey Monday, tell me about cats"
- Wait for green transcript to appear
- Get AI response + voice + YouTube videos

## Features Working
✅ Voice recognition ("Hey Monday" trigger)
✅ ElevenLabs TTS (Sarah voice)  
✅ Perplexity Sonar API responses
✅ YouTube video search with real videos
✅ Three modes: Basic, Reasoning, Deep Research
✅ Clean voice output (no overlapping speech)

## API Keys Configured
- Perplexity API: ✅ Configured
- ElevenLabs API: ✅ Configured  
- YouTube API: ✅ Configured

## Voice Commands
- "Hey Monday, [question]" - Basic mode
- "Hey Monday, think about [topic]" - Reasoning mode  
- "Hey Monday, research into [topic]" - Deep research mode
- "Goodbye" - End conversation

## Browser Requirements
- Chrome/Edge recommended
- Allow microphone permissions
- HTTPS or localhost required for voice

## Troubleshooting
1. Refresh browser (Ctrl+F5)
2. Allow microphone permissions
3. Check both terminals are running
4. Backend must start BEFORE frontend 