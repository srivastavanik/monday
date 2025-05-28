#!/bin/bash

# Monday: Learning, Amalgamated - Setup Script
echo "Setting up Monday: Learning, Amalgamated..."

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if required tools are installed
check_dependencies() {
    echo -e "${YELLOW}Checking dependencies...${NC}"
    
    if ! command -v node &> /dev/null; then
        echo -e "${RED}Node.js is not installed. Please install Node.js 18+ and try again.${NC}"
        exit 1
    fi
    
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}Docker is not installed. Please install Docker and try again.${NC}"
        exit 1
    fi
    
    if ! command -v docker-compose &> /dev/null; then
        echo -e "${RED}Docker Compose is not installed. Please install Docker Compose and try again.${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}All dependencies found!${NC}"
}

# Create SSL certificates for local development (required for WebXR)
create_ssl_certificates() {
    echo -e "${YELLOW}Creating SSL certificates for local development...${NC}"
    
    mkdir -p frontend/certs
    
    # Generate self-signed certificate
    openssl req -x509 -newkey rsa:4096 -keyout frontend/certs/key.pem -out frontend/certs/cert.pem -days 365 -nodes -subj "/C=US/ST=State/L=City/O=Organization/CN=localhost"
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}SSL certificates created successfully!${NC}"
        echo -e "${YELLOW}Note: You'll need to accept the self-signed certificate in your browser.${NC}"
    else
        echo -e "${RED}Failed to create SSL certificates. WebXR may not work properly.${NC}"
    fi
}

# Set up environment variables
setup_environment() {
    echo -e "${YELLOW}Setting up environment variables...${NC}"
    
    if [ ! -f .env ]; then
        cp env.example .env
        echo -e "${GREEN}Environment file created from template.${NC}"
        echo -e "${YELLOW}Please edit .env with your API keys before starting the application.${NC}"
    else
        echo -e "${YELLOW}Environment file already exists.${NC}"
    fi
}

# Install dependencies
install_dependencies() {
    echo -e "${YELLOW}Installing dependencies...${NC}"
    
    # Root dependencies
    npm install
    
    # Frontend dependencies
    echo -e "${YELLOW}Installing frontend dependencies...${NC}"
    cd frontend && npm install && cd ..
    
    # Backend dependencies
    echo -e "${YELLOW}Installing backend dependencies...${NC}"
    cd backend && npm install && cd ..
    
    echo -e "${GREEN}Dependencies installed successfully!${NC}"
}

# Create necessary directories
create_directories() {
    echo -e "${YELLOW}Creating necessary directories...${NC}"
    
    mkdir -p logs
    mkdir -p data/postgres
    mkdir -p data/redis
    mkdir -p data/neo4j
    
    echo -e "${GREEN}Directories created!${NC}"
}

# Main setup function
main() {
    echo -e "${GREEN}Monday: Learning, Amalgamated Setup${NC}"
    echo "=================================="
    
    check_dependencies
    setup_environment
    create_ssl_certificates
    create_directories
    install_dependencies
    
    echo ""
    echo -e "${GREEN}Setup completed successfully!${NC}"
    echo ""
    echo "Next steps:"
    echo "1. Edit .env with your API keys (Perplexity, ElevenLabs)"
    echo "2. Start the application with: npm run docker:up"
    echo "3. Or start in development mode with: npm run dev"
    echo ""
    echo "For Quest development:"
    echo "1. Enable Developer Mode on your Quest headset"
    echo "2. Connect to the same WiFi network"
    echo "3. Navigate to https://YOUR_IP_ADDRESS in Quest browser"
    echo "4. Accept the self-signed certificate warning"
    echo "5. Allow microphone permissions when prompted"
    echo ""
    echo -e "${YELLOW}Happy learning with Monday!${NC}"
}

# Run main function
main 