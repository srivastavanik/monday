# Monday: Learning, Amalgamated - Windows Setup Script
Write-Host "Setting up Monday: Learning, Amalgamated..." -ForegroundColor Green

# Function to check if a command exists
function Test-Command($cmdname) {
    return [bool](Get-Command -Name $cmdname -ErrorAction SilentlyContinue)
}

# Check dependencies
Write-Host "Checking dependencies..." -ForegroundColor Yellow

if (-not (Test-Command node)) {
    Write-Host "Node.js is not installed. Please install Node.js 18+ and try again." -ForegroundColor Red
    exit 1
}

if (-not (Test-Command docker)) {
    Write-Host "Docker is not installed. Please install Docker Desktop and try again." -ForegroundColor Red
    exit 1
}

Write-Host "All dependencies found!" -ForegroundColor Green

# Create SSL certificates for local development
Write-Host "Creating SSL certificates for local development..." -ForegroundColor Yellow

New-Item -ItemType Directory -Force -Path "frontend\certs" | Out-Null

if (Test-Command mkcert) {
    Write-Host "Using mkcert to generate SSL certificates..." -ForegroundColor Cyan
    # Navigate to the certs directory
    Push-Location "frontend\certs"
    
    # Generate certificates
    mkcert localhost 127.0.0.1 ::1
    
    # Check if certificate generation was successful
    if (Test-Path ".\localhost+2.pem") { # mkcert generates files like localhost+2.pem and localhost+2-key.pem
        Rename-Item -Path ".\localhost+2.pem" -NewName "cert.pem"
        Rename-Item -Path ".\localhost+2-key.pem" -NewName "key.pem"
        Write-Host "SSL certificates created successfully using mkcert!" -ForegroundColor Green
        Write-Host "Certificates were automatically added to your local trust store." -ForegroundColor Green
    } else {
        Write-Host "Failed to create SSL certificates using mkcert. WebXR may not work properly." -ForegroundColor Red
    }
    # Return to the previous directory
    Pop-Location
} else {
    Write-Host "mkcert is not installed or not found in PATH." -ForegroundColor Red
    Write-Host "Please install mkcert (https://github.com/FiloSottile/mkcert) and ensure it's in your PATH." -ForegroundColor Red
    Write-Host "Falling back to OpenSSL or placeholder certificates..." -ForegroundColor Yellow
    if (Test-Command openssl) {
        # Generate self-signed certificate using OpenSSL
        openssl req -x509 -newkey rsa:4096 -keyout frontend\certs\key.pem -out frontend\certs\cert.pem -days 365 -nodes -subj "/C=US/ST=State/L=City/O=Organization/CN=localhost"
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "SSL certificates created successfully using OpenSSL!" -ForegroundColor Green
            Write-Host "Note: You'll need to accept the self-signed certificate in your browser." -ForegroundColor Yellow
        } else {
            Write-Host "Failed to create SSL certificates using OpenSSL. WebXR may not work properly." -ForegroundColor Red
        }
    } else {
        Write-Host "OpenSSL not found. Creating placeholder certificates..." -ForegroundColor Yellow
        "# Placeholder certificate - install mkcert or OpenSSL to generate real certs" | Out-File -FilePath "frontend\certs\cert.pem"
        "# Placeholder key - install mkcert or OpenSSL to generate real certs" | Out-File -FilePath "frontend\certs\key.pem"
    }
}

# Set up environment variables
Write-Host "Setting up environment variables..." -ForegroundColor Yellow

if (-not (Test-Path ".env")) {
    Copy-Item "env.example" ".env"
    Write-Host "Environment file created from template." -ForegroundColor Green
    Write-Host "Please edit .env with your API keys before starting the application." -ForegroundColor Yellow
} else {
    Write-Host "Environment file already exists." -ForegroundColor Yellow
}

# Install dependencies
Write-Host "Installing dependencies..." -ForegroundColor Yellow

# Root dependencies
npm install

if ($LASTEXITCODE -ne 0) {
    Write-Host "Failed to install root dependencies" -ForegroundColor Red
    exit 1
}

# Frontend dependencies
Write-Host "Installing frontend dependencies..." -ForegroundColor Yellow
Set-Location frontend
npm install
if ($LASTEXITCODE -ne 0) {
    Write-Host "Failed to install frontend dependencies" -ForegroundColor Red
    exit 1
}
Set-Location ..

# Backend dependencies
Write-Host "Installing backend dependencies..." -ForegroundColor Yellow
Set-Location backend
npm install
if ($LASTEXITCODE -ne 0) {
    Write-Host "Failed to install backend dependencies" -ForegroundColor Red
    exit 1
}
Set-Location ..

Write-Host "Dependencies installed successfully!" -ForegroundColor Green

# Create necessary directories
Write-Host "Creating necessary directories..." -ForegroundColor Yellow

New-Item -ItemType Directory -Force -Path "logs" | Out-Null
New-Item -ItemType Directory -Force -Path "data\postgres" | Out-Null
New-Item -ItemType Directory -Force -Path "data\redis" | Out-Null
New-Item -ItemType Directory -Force -Path "data\neo4j" | Out-Null

Write-Host "Directories created!" -ForegroundColor Green

Write-Host ""
Write-Host "Setup completed successfully!" -ForegroundColor Green
Write-Host ""
Write-Host "Next steps:"
Write-Host "1. Edit .env with your API keys (Perplexity, ElevenLabs)"
Write-Host "2. Start the application with: npm run docker:up"
Write-Host "3. Or start in development mode with: npm run dev"
Write-Host ""
Write-Host "For Quest development:"
Write-Host "1. Enable Developer Mode on your Quest headset"
Write-Host "2. Connect to the same WiFi network"
Write-Host "3. Find your IP address with: ipconfig"
Write-Host "4. Navigate to https://YOUR_IP_ADDRESS in Quest browser"
Write-Host "5. Accept the self-signed certificate warning"
Write-Host "6. Allow microphone permissions when prompted"
Write-Host ""
Write-Host "Happy learning with Monday!" -ForegroundColor Yellow 