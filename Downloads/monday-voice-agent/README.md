# Monday Voice Agent

A Next.js-based voice-enabled AI assistant interface for the Monday Learning System.

## Features

- Voice recognition and text-to-speech capabilities
- Integration with Perplexity AI models (sonar-pro, sonar-reasoning-pro, sonar-deep-research)
- 3D spatial interface with interactive panels
- Real-time WebSocket communication with backend
- Progressive reasoning and research display
- Intelligent voice feedback loop prevention

## Getting Started

1. Install dependencies:
```bash
pnpm install
```

2. Run the development server:
```bash
pnpm dev
```

3. Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

## Tech Stack

- Next.js 15
- React 19
- TypeScript
- Tailwind CSS
- Three.js / React Three Fiber
- Radix UI components
- WebSocket for real-time communication

## Voice System

The voice system includes:
- Speech recognition with intelligent filtering
- Text-to-speech with audio context management
- Microphone control with feedback loop prevention
- Command validation and voice signature filtering

## Backend Integration

Connects to the Monday Learning System backend for:
- AI model selection and routing
- Progressive reasoning and research
- Real-time response streaming
- Voice command processing 