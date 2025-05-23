# Monday: Learning, Amalgamated

A voice-driven, ambient AI learning companion that operates in virtual reality, specifically designed for Oculus Quest 2/3 headsets via WebXR. The application leverages Perplexity's Sonar API ecosystem to create an immersive educational environment where knowledge materializes in three-dimensional space around the user.

## ✅ Complete Implementation Status

**All components are fully built and functional:**

### Frontend (React + WebXR)
- ✅ Complete React/TypeScript application with WebXR support
- ✅ Zustand store for state management
- ✅ Voice recognition with Web Speech API  
- ✅ WebSocket real-time communication
- ✅ 3D information panels with Perplexity branding
- ✅ Spatial orchestration and layout management
- ✅ Performance monitoring for Quest optimization
- ✅ HTTPS configuration for WebXR requirements

### Backend (Node.js API)
- ✅ Express server with Socket.IO integration
- ✅ Complete Perplexity Sonar API integration (Basic, Reasoning, Deep Research)
- ✅ Comprehensive voice command parsing
- ✅ Session management with Redis
- ✅ Database connections (PostgreSQL, Redis, Neo4j)
- ✅ Analytics and performance tracking
- ✅ Security middleware and rate limiting

### Python Services
- ✅ FastAPI microservices for ML processing
- ✅ Research synthesis endpoints
- ✅ Knowledge graph processing
- ✅ Visualization generation
- ✅ Learning path optimization

### Infrastructure
- ✅ Docker containerization for all services
- ✅ Nginx reverse proxy with WebXR optimization
- ✅ SSL certificate generation for local development
- ✅ Database initialization and table creation

## Quick Start

### Prerequisites
- Node.js 18+
- Docker and Docker Compose
- Oculus Quest 2/3 headset (for VR testing)
- **Required API Keys:**
  - Perplexity API key (get from [Perplexity AI](https://www.perplexity.ai/))
  - ElevenLabs API key (get from [ElevenLabs](https://elevenlabs.io/))

### Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd monday-learning-amalgamated
   ```

2. **Run setup script**
   
   **On Windows:**
   ```powershell
   .\scripts\setup.ps1
   ```
   
   **On macOS/Linux:**
   ```bash
   chmod +x scripts/setup.sh
   ./scripts/setup.sh
   ```

3. **🔑 Configure your API keys**
   
   Edit the `.env` file in the root directory:
   ```bash
   # Required API keys
   PERPLEXITY_API_KEY=your_perplexity_api_key_here
   ELEVENLABS_API_KEY=your_elevenlabs_api_key_here
   
   # Optional: Customize other settings
   ELEVENLABS_VOICE_ID=21m00Tcm4TlvDq8ikWAM
   NODE_ENV=development
   # ... other configuration options
   ```

4. **Start the application**
   
   **With Docker (Recommended):**
   ```bash
   npm run docker:up
   ```
   
   **Development mode:**
   ```bash
   npm run dev
   ```

5. **Access the application**
   - **Desktop**: Navigate to `https://localhost:3000` (accept SSL warning)
   - **Quest VR**: Navigate to `https://YOUR_IP_ADDRESS` in Quest browser

### VR Setup for Quest

1. **Enable Developer Mode** on your Quest headset
2. **Connect to same WiFi** network as your development machine
3. **Find your IP address:**
   - Windows: `ipconfig`
   - macOS/Linux: `ifconfig`
4. **Navigate to `https://YOUR_IP_ADDRESS`** in Quest browser
5. **Accept SSL certificate warning** (required for WebXR)
6. **Allow microphone permissions** when prompted
7. **Say "Monday, hello"** to begin

## Features

### Voice Commands
- **Basic Query**: "Monday, [question]"
- **Reasoning Mode**: "Monday, think about [topic]"
- **Deep Research**: "Monday, deep dive on [topic]"
- **Spatial Control**: "Monday, bring this closer" / "Monday, push that away"
- **Focus Mode**: "Monday, focus mode"
- **Context Reference**: "Monday, what about [topic] from earlier?"

### Learning Modes
1. **Basic**: Quick answers using Perplexity Sonar Basic API
2. **Reasoning**: Step-by-step analysis with visible thought processes
3. **Deep Research**: Multi-source research with expanding knowledge webs

### 3D Visualizations
- Information panels arrange automatically in VR space
- Reasoning chains visualize logical connections
- Research webs show source relationships
- Focus mode dims distractions for concentration

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

### Infrastructure
- **Containerization**: Docker & Docker Compose
- **Reverse Proxy**: Nginx with WebXR optimization
- **SSL**: Self-signed certificates for local HTTPS

## Development

### Project Structure
```
monday-learning-amalgamated/
├── .env                    # 🔑 API keys and configuration
├── frontend/              # WebXR React application
│   ├── src/
│   │   ├── components/    # VR components and UI
│   │   ├── hooks/         # Voice recognition, WebSocket
│   │   ├── store/         # Zustand state management
│   │   └── types/         # TypeScript definitions
│   └── certs/             # SSL certificates for HTTPS
├── backend/               # Node.js API server
│   ├── src/
│   │   ├── routes/        # API endpoints
│   │   ├── services/      # Perplexity integration
│   │   ├── sockets/       # WebSocket handlers
│   │   └── database/      # Database connections
├── python-services/       # ML and research synthesis
├── nginx/                 # Reverse proxy configuration
├── scripts/               # Setup and utility scripts
└── docker-compose.yml     # Container orchestration
```

### Environment Variables

The `.env` file contains all configuration:

```bash
# API Keys (Required)
PERPLEXITY_API_KEY=your_key_here
ELEVENLABS_API_KEY=your_key_here

# Database URLs
DATABASE_URL=postgresql://monday_user:monday_password@localhost:5432/monday_learning
REDIS_URL=redis://localhost:6379
NEO4J_URI=bolt://localhost:7687

# Security (Change in production)
JWT_SECRET=monday_jwt_secret_change_in_production
SESSION_SECRET=monday_session_secret_change_in_production

# Performance Settings
MAX_CONCURRENT_PANELS=8
TARGET_FPS=72
```

### Performance Targets

- **72 FPS** minimum on Quest hardware
- **Sub-200ms** response time for voice queries
- **Max 8 concurrent** information panels
- **150k triangle limit** in VR view
- **Automatic quality scaling** based on performance

## API Integration

### Perplexity Sonar
- **Basic** (`llama-3.1-sonar-small-128k-online`): Quick queries
- **Reasoning** (`llama-3.1-sonar-large-128k-online`): Step-by-step analysis
- **Deep Research** (`llama-3.1-sonar-huge-128k-online`): Multi-source research

### ElevenLabs TTS
- Natural voice synthesis for Monday's responses
- Optimized for clarity and VR audio

## Troubleshooting

### Common Issues

1. **"Cannot find .env file"**
   - Run the setup script: `.\scripts\setup.ps1` (Windows) or `./scripts/setup.sh` (Linux/Mac)
   - The setup script creates `.env` from `env.example`

2. **WebXR not working**
   - Ensure HTTPS is enabled (required for WebXR)
   - Accept self-signed certificate warning
   - Check Quest developer mode is enabled

3. **Voice recognition failing**
   - Check microphone permissions in browser
   - Ensure HTTPS connection (required for microphone access)
   - Try different browsers (Chrome recommended)

4. **API errors**
   - Verify API keys are correctly set in `.env`
   - Check internet connectivity
   - Review logs: `docker-compose logs backend`

### Getting Help

1. Check the development guide: `DEVELOPMENT.md`
2. Review the demo script: `DEMO.md`
3. Check GitHub issues
4. Ensure all prerequisites are installed

## Contributing

1. Fork the repository
2. Create a feature branch
3. Follow TypeScript strict mode
4. Add tests for new features
5. Test on actual Quest hardware

## License

MIT License - see LICENSE file for details

---

**Ready to learn with Monday? Run the setup script and start exploring knowledge in VR!** 🚀 