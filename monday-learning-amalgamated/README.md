# Monday: Learning, Amalgamated

A voice-driven, ambient AI learning companion that operates in virtual reality, specifically designed for Oculus Quest 2/3 headsets via WebXR. The application leverages Perplexity's Sonar API ecosystem to create an immersive educational environment where knowledge materializes in three-dimensional space around the user.

## Features

- **Voice-First Interaction**: Operates primarily through voice commands with "Monday" activation
- **Spatial Information Management**: Information panels arrange automatically in 3D space
- **Three Learning Modes**:
  - Basic queries using Perplexity Sonar API
  - Reasoning mode with visible thought processes
  - Deep research with expanding knowledge webs
- **3D Visualizations**: Pre-built models for computer science concepts and dynamic generation
- **Perplexity Brand Integration**: Consistent visual theme and typography
- **AI Personality**: Professional yet approachable assistant with emotional intelligence

## Technology Stack

### Frontend
- **WebXR**: A-Frame framework with Three.js
- **React**: With react-xr for UI overlay components
- **Voice**: Web Speech API for recognition
- **Audio**: Web Audio API for spatial audio
- **Communication**: WebSocket for real-time updates

### Backend
- **Runtime**: Node.js with Express
- **Database**: PostgreSQL, Redis, Neo4j
- **APIs**: Perplexity Sonar, ElevenLabs TTS
- **Services**: Python microservices for ML processing

## Quick Start

### Prerequisites
- Node.js 18+
- Docker and Docker Compose
- Oculus Quest 2/3 headset
- Perplexity API key
- ElevenLabs API key

### Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd monday-learning-amalgamated
   ```

2. **Set up environment variables**
   ```bash
   cp env.example .env
   # Edit .env with your API keys
   ```

3. **Install dependencies**
   ```bash
   npm run install-all
   ```

4. **Start with Docker**
   ```bash
   npm run docker:up
   ```

5. **Or start development servers**
   ```bash
   npm run dev
   ```

### VR Setup

1. Enable Developer Mode on your Quest headset
2. Connect to the same WiFi network as your development machine
3. Navigate to your machine's IP address in the Quest browser
4. Allow microphone permissions when prompted
5. Say "Monday, hello" to begin

## Usage

### Voice Commands

- **Basic Query**: "Monday, [question]"
- **Reasoning Mode**: "Monday, think about [topic]"
- **Deep Research**: "Monday, deep dive on [topic]"
- **Spatial Control**: "Monday, bring this closer" / "Monday, push that away"
- **Focus Mode**: "Monday, focus mode"
- **Context Reference**: "Monday, what about [topic] from earlier?"

### Spatial Navigation

- Look around to see information arranged in 180-degree field of view
- Active content appears in center, related information on sides
- Voice commands control spatial arrangement automatically
- Focus mode dims non-active panels for concentration

## Development

### Project Structure
```
monday-learning-amalgamated/
├── frontend/          # WebXR React application
├── backend/           # Node.js API server
├── python-services/   # ML and research synthesis services
├── nginx/            # Reverse proxy configuration
└── docker-compose.yml # Container orchestration
```

### API Integration

The application integrates with:
- **Perplexity Sonar API**: Basic, Reasoning, and Deep Research modes
- **ElevenLabs API**: Text-to-speech with personality
- **WebXR Device API**: Quest headset integration

### Performance Targets

- Maintain 72 FPS on Quest hardware
- Support up to 8 concurrent information panels
- Maximum 150,000 triangles in view
- Sub-200ms response time for voice queries

## Deployment

### Production Build
```bash
npm run build
npm run docker:build
```

### Environment Setup
- Configure reverse proxy for HTTPS
- Set up SSL certificates for WebXR
- Configure API rate limiting
- Enable production logging

## Contributing

1. Fork the repository
2. Create a feature branch
3. Follow the existing code style
4. Add tests for new features
5. Submit a pull request

## License

MIT License - see LICENSE file for details

## Support

For issues and questions:
- Check the GitHub issues
- Review the troubleshooting guide
- Contact the development team 