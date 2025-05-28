# Monday: 

A voice-driven AI learning companion that operates in virtual reality, specifically designed for Oculus Quest 2/3 headsets via WebXR. The application leverages Perplexity's Sonar API ecosystem to create an immersive educational environment where knowledge materializes around the user.

## Implementation

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
  - Youtube API key (get from [YouTube Data API](https://developers.google.com/youtube/v3/getting-started/))

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
7. **Say "Hey Monday"** to begin

## Features

### Voice Commands
- **Basic Query**: "Monday, [question]"
- **Reasoning Mode**: "Monday, think about [topic]"
- **Deep Research**: "Monday, research into [topic]"

### Learning Modes
1. **Basic**: Quick answers using Perplexity Sonar Basic API
2. **Reasoning**: Step-by-step analysis with visible thought processes
3. **Deep Research**: Multi-source research with expanding knowledge panels

### 3D Visualizations
- Information panels arrange automatically in VR space
- Reasoning chains visualize logical connections
- Research webs show source relationships
- Focus mode dims distractions for concentration

## Technology Stack

### Frontend  
- **React + TypeScript**: Modern, type-safe UI with Tailwind CSS styling.
- **WebXR + Three.js + A-Frame**: For immersive 3D visualizations and optionally VR-enabled experiences.
- **React-XR**: UI overlay components in 3D space.
- **Voice Recognition**: Web Speech API powers hands-free, voice-first interaction.
- **Spatial Audio**: Web Audio API creates directional sound tied to virtual panels.
- **Real-Time Communication**: WebSocket-based connection for low-latency interaction and real-time updates.

### Backend

- **Node.js + Express**  
  REST and WebSocket server managing requests and streaming updates.

- **Perplexity Sonar API**  
  AI-generated answers tailored to different learning modes (Basic, Reasoning, Deep Research).

- **ElevenLabs API**  
  Natural-sounding text-to-speech responses.

- **Database Layer**  
  - **Neo4j**: For reasoning mode (step-by-step logic trees).


### Infrastructure
- **Reverse Proxy**: Nginx with WebXR optimization
- **SSL**: Self-signed certificates for local HTTPS

## Development

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

2. **Voice recognition failing**
   - Check microphone permissions in browser
   - Ensure HTTPS connection (required for microphone access)
   - Try different browsers (Chrome recommended)

3. **API errors**
   - Verify API keys are correctly set in `.env`
   - Check internet connectivity
   - Review logs: `docker-compose logs backend`


## Contributing

1. Fork the repository
2. Create a feature branch
3. Follow TypeScript strict mode
4. Add tests for new features
5. Test on actual Quest hardware


**Ready to learn with Monday? Run the setup script and start exploring knowledge!** 
