# Monday: Learning, Amalgamated - Development Guide

## Quick Start

### Prerequisites
- Node.js 18+
- Docker & Docker Compose
- Oculus Quest 2/3 headset (for testing)
- Perplexity API key
- ElevenLabs API key

### Setup

1. **Clone and setup the project:**
   ```bash
   git clone <repository-url>
   cd monday-learning-amalgamated
   chmod +x scripts/setup.sh
   ./scripts/setup.sh
   ```

2. **Configure environment variables:**
   ```bash
   cp env.example .env
   # Edit .env with your API keys
   ```

3. **Start the application:**
   ```bash
   # With Docker (recommended)
   npm run docker:up
   
   # Or development mode
   npm run dev
   ```

4. **Access the application:**
   - Desktop: `https://localhost:3000` (accept SSL warning)
   - Quest: `https://YOUR_IP_ADDRESS` (enable developer mode first)

## Architecture Overview

### Frontend (React + WebXR)
- **Location**: `frontend/`
- **Framework**: React 18 with TypeScript
- **3D Engine**: Three.js + A-Frame + @react-three/fiber
- **VR**: WebXR API with @react-three/xr
- **Voice**: Web Speech API
- **State**: Zustand store
- **Build**: Vite with PWA support

### Backend (Node.js API)
- **Location**: `backend/`
- **Runtime**: Node.js 18+ with TypeScript
- **Framework**: Express.js
- **Real-time**: Socket.IO for WebSocket communication
- **APIs**: Perplexity Sonar integration
- **Databases**: PostgreSQL, Redis, Neo4j

### Python Services
- **Location**: `python-services/`
- **Framework**: FastAPI
- **Purpose**: ML processing, research synthesis
- **AI/ML**: Natural language processing, knowledge graphs

### Infrastructure
- **Reverse Proxy**: Nginx
- **Containerization**: Docker & Docker Compose
- **SSL**: Self-signed certificates for local HTTPS (required for WebXR)

## Development Workflow

### Frontend Development

1. **Start development server:**
   ```bash
   cd frontend
   npm run dev
   ```

2. **Key files:**
   - `src/App.tsx` - Main application component
   - `src/components/MondayScene.tsx` - 3D VR scene
   - `src/hooks/useVoiceRecognition.ts` - Voice command handling
   - `src/store/mondayStore.ts` - Global state management

3. **VR Testing:**
   - Enable Quest Developer Mode
   - Connect to same WiFi as development machine
   - Navigate to `https://YOUR_IP_ADDRESS:3000`
   - Accept SSL certificate warning
   - Allow microphone permissions

### Backend Development

1. **Start development server:**
   ```bash
   cd backend
   npm run dev
   ```

2. **Key files:**
   - `src/index.ts` - Main server entry point
   - `src/services/perplexity.ts` - Perplexity API integration
   - `src/sockets/voiceHandler.ts` - Voice command processing
   - `src/routes/` - API endpoints

3. **API Testing:**
   - Health check: `GET http://localhost:3001/health`
   - API docs: Available via OpenAPI/Swagger (add if needed)

### Python Services Development

1. **Start development server:**
   ```bash
   cd python-services
   pip install -r requirements.txt
   uvicorn main:app --reload --host 0.0.0.0 --port 8000
   ```

2. **API Documentation:**
   - FastAPI docs: `http://localhost:8000/docs`
   - ReDoc: `http://localhost:8000/redoc`

## Voice Commands

Monday responds to these voice patterns:

- **Basic Query**: "Monday, [question]"
- **Reasoning Mode**: "Monday, think about [topic]"
- **Deep Research**: "Monday, deep dive on [topic]"
- **Spatial Control**: "Monday, bring this closer" / "Monday, push that away"
- **Focus Mode**: "Monday, focus mode"
- **Context Reference**: "Monday, what about [topic] from earlier?"

## WebXR Development Tips

### Quest-Specific Considerations

1. **Performance Targets:**
   - Maintain 72 FPS minimum
   - Keep draw calls under 100
   - Limit triangles to 150,000 in view
   - Monitor memory usage (target <500MB)

2. **Input Methods:**
   - Primary: Voice commands
   - Secondary: Hand tracking
   - Fallback: Controller input

3. **UI Design:**
   - Use high contrast (Perplexity brand colors)
   - Minimum font size: 16px equivalent in VR
   - Touch targets: minimum 44px equivalent
   - Readable distance: 1-3 meters

### WebXR API Usage

```javascript
// Check VR support
if (navigator.xr) {
  const isSupported = await navigator.xr.isSessionSupported('immersive-vr')
}

// Start VR session
const session = await navigator.xr.requestSession('immersive-vr')
```

### Performance Monitoring

The app includes built-in performance monitoring:
- FPS counter
- Frame time measurement
- Memory usage tracking
- Automatic quality adjustment

## API Integration

### Perplexity Sonar API

Three tiers implemented:
1. **Basic** (`llama-3.1-sonar-small-128k-online`) - Quick queries
2. **Reasoning** (`llama-3.1-sonar-large-128k-online`) - Step-by-step analysis  
3. **Deep Research** (`llama-3.1-sonar-huge-128k-online`) - Multi-source research

### ElevenLabs TTS

Voice synthesis for Monday's responses:
- Voice ID: Configurable in environment
- Settings optimized for clarity and naturalness

## Database Schema

### PostgreSQL
- User sessions and analytics
- Learning progress tracking
- Performance metrics

### Redis
- Session state caching
- Real-time data
- WebSocket session management

### Neo4j
- Knowledge graph storage
- Concept relationships
- Learning path optimization

## Deployment

### Development
```bash
npm run dev  # Start all services in development mode
```

### Production
```bash
npm run build      # Build all services
npm run docker:build  # Build Docker images
npm run docker:up  # Start with Docker Compose
```

### Environment Variables

Required for production:
- `PERPLEXITY_API_KEY` - Perplexity Sonar API access
- `ELEVENLABS_API_KEY` - Text-to-speech service
- `JWT_SECRET` - Session security
- `DATABASE_URL` - PostgreSQL connection
- `REDIS_URL` - Redis connection
- `NEO4J_URI` - Neo4j connection

## Troubleshooting

### Common Issues

1. **WebXR not working:**
   - Ensure HTTPS is enabled (required for WebXR)
   - Check browser WebXR support
   - Verify Quest developer mode is enabled

2. **Voice recognition failing:**
   - Check microphone permissions
   - Ensure HTTPS connection
   - Test with different browsers

3. **Performance issues:**
   - Check FPS counter in development mode
   - Reduce quality settings
   - Monitor memory usage

4. **API errors:**
   - Verify API keys in environment
   - Check network connectivity
   - Review logs for detailed errors

### Debugging

1. **Frontend debugging:**
   ```bash
   # Enable verbose logging
   localStorage.setItem('debug', 'monday:*')
   ```

2. **Backend debugging:**
   ```bash
   # Set log level
   export LOG_LEVEL=debug
   npm run dev
   ```

3. **Quest debugging:**
   - Use Chrome DevTools with Quest connected via USB
   - Enable WebXR debugging in chrome://flags
   - Check Quest browser console logs

## Contributing

1. **Code Style:**
   - TypeScript strict mode
   - ESLint configuration provided
   - Prettier for formatting

2. **Testing:**
   - Jest for unit tests
   - Cypress for E2E tests (add if needed)
   - Manual VR testing required

3. **Performance:**
   - Monitor FPS in VR
   - Profile memory usage
   - Test on actual Quest hardware

## Resources

- [WebXR Specification](https://www.w3.org/TR/webxr/)
- [Perplexity API Docs](https://docs.perplexity.ai/)
- [Three.js Documentation](https://threejs.org/docs/)
- [Quest Developer Hub](https://developer.oculus.com/)
- [A-Frame School](https://aframe.io/aframe-school/)

## Support

For issues and questions:
- Check GitHub Issues
- Review this development guide
- Test with the demo flow
- Ensure all prerequisites are met 