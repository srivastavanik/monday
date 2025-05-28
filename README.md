# Monday - AI Voice Companion

> **Contest Submission Categories: Most Creative/Fun & Best Deep Research**

An intelligent voice-activated AI companion that combines conversational AI, real-time voice synthesis, and dynamic content discovery. Monday responds to natural voice commands, provides thoughtful responses through multiple AI modes, and enhances learning with relevant video content.

## 🎯 **Quick Start**

### Prerequisites
- Node.js 18+ 
- Chrome/Edge browser (for voice recognition)
- Microphone access

### Installation & Setup

1. **Clone and Navigate**
```bash
git clone https://github.com/srivastavanik/monday.git
cd monday
git checkout final
cd nidsmonday
```

2. **Install Dependencies**
```bash
npm install
```

3. **Start Backend Server** (Terminal 1)
```bash
node backend-server.js
```
✅ You should see: `"Monday backend server running on port 3001"`

4. **Start Frontend** (Terminal 2)  
```bash
npm run dev
```
✅ Go to: **http://localhost:3000**

5. **Grant Permissions**
- Allow microphone access when prompted
- Wait for "🎤 Listening for 'Hey Monday'..." message

6. **Start Talking!**
```
Say: "Hey Monday, tell me about quantum physics"
```

## 🎤 **Voice Commands**

| Command Pattern | AI Mode | Example |
|----------------|---------|---------|
| `"Hey Monday, [question]"` | **Basic** | "Hey Monday, what is photosynthesis?" |
| `"Hey Monday, think about [topic]"` | **Reasoning** | "Hey Monday, think about why the sky is blue" |
| `"Hey Monday, research into [topic]"` | **Deep Research** | "Hey Monday, research into renewable energy" |
| `"Goodbye"` | **Exit** | Ends conversation mode |

## 🧠 **AI Modes & Perplexity Integration**

Monday uses **Perplexity's Sonar API** with three distinct modes:

### 1. **Basic Mode** (sonar-small-online)
- **Purpose**: Quick answers and factual information
- **Trigger**: Natural questions without reasoning keywords
- **API Usage**: Fast responses with web search integration
- **Example**: "What is the capital of France?"

### 2. **Reasoning Mode** (sonar-medium-online)  
- **Purpose**: Analytical thinking and step-by-step explanations
- **Trigger**: Keywords like "think about", "analyze", "explain why"
- **API Usage**: Enhanced reasoning capabilities with citations
- **Example**: "Think about why gravity works the way it does"

### 3. **Deep Research Mode** (sonar-large-online)
- **Purpose**: Comprehensive investigation and detailed analysis  
- **Trigger**: Keywords like "research into", "investigate", "deep dive"
- **API Usage**: Extensive research with multiple sources and citations
- **Example**: "Research into the latest developments in quantum computing"

### **Perplexity API Implementation**
- **Real-time web search**: Accesses current information beyond training data
- **Source citations**: Provides credible references for all claims
- **Multi-modal responses**: Combines text, reasoning, and source validation
- **Adaptive complexity**: Automatically adjusts response depth based on query type

## 🔊 **Voice Technology Stack**

### **Speech Recognition**
- **Library**: react-speech-recognition with Web Speech API
- **Features**: Continuous listening, wake word detection ("Hey Monday")
- **Browser Support**: Chrome, Edge, Safari (with permissions)

### **Text-to-Speech**
- **Service**: ElevenLabs API with "Sarah" voice model
- **Features**: Natural-sounding speech synthesis
- **Voice ID**: `EXAVITQu4vr4xnSDxMaL` (Sarah - clear, professional tone)
- **Optimization**: Single voice stream to prevent overlapping audio

## 🎥 **Dynamic Video Discovery**

### **YouTube Integration**
- **API**: YouTube Data API v3 with intelligent search
- **Keyword Extraction**: Analyzes AI responses to find relevant topics
- **Video Validation**: Ensures all videos are embeddable and accessible
- **Display**: 1 main video player + 2 related suggestions

### **Smart Search Algorithm**
1. **Extract keywords** from AI response content using NLP
2. **Build enhanced queries** based on detected AI mode
3. **Filter results** for educational, embeddable content  
4. **Validate availability** before display

**Example Flow:**
```
User: "Hey Monday, think about photosynthesis"
AI Response: "Photosynthesis involves chlorophyll, sunlight, and carbon dioxide..."
Keywords Extracted: ["photosynthesis", "chlorophyll", "sunlight"]
YouTube Query: "photosynthesis chlorophyll sunlight explained tutorial analysis"
Result: 3 relevant educational videos about photosynthesis
```

## 🏗️ **Architecture**

### **Frontend** (`nidsmonday/`)
- **Framework**: Next.js 15 with TypeScript
- **Styling**: Tailwind CSS with custom gradients
- **Components**: Modular React components for voice, visualization, and video
- **State Management**: React hooks with WebSocket integration

### **Backend** (`nidsmonday/backend-server.js`)
- **Runtime**: Node.js with Express
- **WebSockets**: Socket.IO for real-time communication
- **APIs**: Perplexity, ElevenLabs, YouTube Data API
- **Features**: Mode detection, response processing, error handling

### **Real-time Communication**
```
Browser ←→ WebSocket ←→ Node.js Backend ←→ External APIs
   ↓           ↓              ↓              ↓
Voice Input → Processing → AI Response → Voice Output
   ↓           ↓              ↓              ↓  
Transcript → Mode Detection → Content → Video Search
```

## 🎨 **UI Components**

### **1. Adaptive Visualization Panel**
- **Purpose**: Visual representation of AI thinking process
- **Features**: Dynamic 3D visualizations, progress tracking
- **Mode Awareness**: Different visualizations per AI mode

### **2. Voice Processing Panel** 
- **Purpose**: Central conversation interface
- **Features**: Live transcript, response display, mode indicators
- **Real-time**: Shows typing effect and voice status

### **3. YouTube Discovery Panel**
- **Purpose**: Educational video enhancement
- **Features**: Smart video search, embedded player, related content
- **Integration**: Keyword-driven search based on AI responses

## 🔧 **Configuration**

All API keys are pre-configured in the codebase:

- **Perplexity API**: `pplx-CwPQD...` (Sonar models)
- **ElevenLabs**: `sk_6be41ddad8e1f62b...` (Sarah voice)  
- **YouTube Data API**: `AIzaSyDn_zCV8AGj...` (Video search)

## 🏆 **Contest Categories**

### **Most Creative/Fun**
- **Natural Voice Interaction**: Conversational AI with wake words
- **Multi-Modal Experience**: Voice + Visual + Video integration
- **Personality**: "Monday" character with consistent voice and responses
- **User Experience**: Seamless voice-to-AI-to-content pipeline

### **Best Deep Research**
- **Perplexity Integration**: Advanced AI models for research
- **Source Validation**: Real citations and credible references
- **Adaptive Complexity**: Different research depths (Basic/Reasoning/Deep)
- **Knowledge Enhancement**: AI responses augmented with educational videos

## 🚀 **Features**

✅ **Voice Recognition** - "Hey Monday" wake word detection  
✅ **Multi-Mode AI** - Basic, Reasoning, Deep Research modes  
✅ **Natural Speech** - ElevenLabs voice synthesis  
✅ **Smart Video Search** - AI-driven YouTube content discovery  
✅ **Real-time Processing** - WebSocket-based communication  
✅ **Responsive UI** - Modern, accessible interface  
✅ **Error Handling** - Robust fallbacks and user feedback  
✅ **Citation Support** - Source references for all AI responses  

## 🐛 **Troubleshooting**

| Issue | Solution |
|-------|----------|
| No voice recognition | Refresh browser (Ctrl+F5), allow microphone |
| "Backend not connected" | Ensure `node backend-server.js` is running first |
| No voice output | Check ElevenLabs API key and audio permissions |
| Videos not loading | Verify YouTube API key and internet connection |
| Mode not switching | Check console for transcript processing logs |

This is a hackathon submission to Perplexity Sonar (Devpost). Aiming for Most Creative / Fun + Best Deep Research.

---

**Built with ❤️ for the AI Voice Assistant Contest**  
*Demonstrating the future of conversational AI and multimodal learning* 
