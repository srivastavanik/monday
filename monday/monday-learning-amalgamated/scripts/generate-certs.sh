#!/bin/bash

# Generate SSL certificates for Monday WebXR development
echo "Generating SSL certificates for Monday WebXR development..."

# Create certs directory
mkdir -p frontend/certs

# Generate private key
openssl genrsa -out frontend/certs/key.pem 4096

# Generate certificate signing request
openssl req -new -key frontend/certs/key.pem -out frontend/certs/csr.pem -subj "/C=US/ST=CA/L=SF/O=Monday/CN=localhost"

# Generate self-signed certificate
openssl x509 -req -days 365 -in frontend/certs/csr.pem -signkey frontend/certs/key.pem -out frontend/certs/cert.pem

# Clean up CSR
rm frontend/certs/csr.pem

# Set appropriate permissions
chmod 600 frontend/certs/key.pem
chmod 644 frontend/certs/cert.pem

echo "SSL certificates generated successfully!"
echo "Location: frontend/certs/"
echo "- Certificate: cert.pem"
echo "- Private Key: key.pem"
echo ""
echo "Note: You'll need to accept the self-signed certificate warning in your browser."
echo "This is required for WebXR to work with microphone access." 