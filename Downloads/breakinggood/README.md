# Breaking Good

A comprehensive neuropharmacology drug discovery platform for designing, simulating, and analyzing novel ADHD therapeutics.

## Overview

Breaking Good is an end-to-end drug discovery platform that integrates advanced computational chemistry tools, AI-powered molecular design, and collaborative research capabilities. The platform focuses on designing novel therapeutics for ADHD with improved efficacy and reduced side effects compared to existing treatments.

## Features

- **AI-Powered Molecule Design**: Leverages Claude by Anthropic to generate and refine novel molecular structures with detailed reasoning
- **Real-time 3D Molecule Visualization**: Interactive visualization of molecular structures using 3Dmol.js
- **Advanced Simulations**: Molecular docking, dynamics, and ADMET prediction capabilities
- **Database Integration**: Search and query PubMed and ChEMBL databases
- **Collaboration Tools**: Team-based project management with real-time commenting and activity tracking
- **Regulatory Analysis**: Pharmaceutical development timeline prediction and regulatory pathway analysis

## Prerequisites

- Node.js (v14+)
- Python (3.8+) with RDKit installed
- Anthropic API key for Claude integration

## Installation

### Backend Setup

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/breaking-good.git
   cd breaking-good/backend
   ```

2. Install dependencies:
   ```
   npm install
   ```

3. Create a `.env` file in the backend directory with the following variables:
   ```
   PORT=5000
   JWT_SECRET=your_secret_key_for_jwt
   ANTHROPIC_API_KEY=your_anthropic_api_key
   ```

4. Install Python dependencies:
   ```
   pip install -r requirements.txt
   ```

### Frontend Setup

1. Navigate to the frontend directory:
   ```
   cd ../frontend
   ```

2. Install dependencies:
   ```
   npm install
   ```

## Running the Platform

### Development Mode

1. Start the backend server:
   ```
   cd backend
   npm run dev
   ```

2. In a separate terminal, start the frontend development server:
   ```
   cd frontend
   npm start
   ```

3. Open your browser and navigate to `http://localhost:3000`

### Production Mode

1. Build the frontend:
   ```
   cd frontend
   npm run build
   ```

2. Start the backend server:
   ```
   cd backend
   npm start
   ```

3. Access the application at `http://localhost:5000`

## User Guide

### Getting Started

1. Create an account or log in with the default admin credentials:
   - Username: admin
   - Password: admin123

2. Navigate to the Molecule Designer to start creating new molecular designs

### Using AI-Powered Molecule Design

1. Navigate to the "AI Designer" tab
2. Enter your design requirements, specifying target receptors and constraints
3. Click "Generate with AI" to have Claude create novel molecular designs
4. Explore the detailed thinking process to understand the reasoning behind each design
5. Save promising molecules for further analysis

### Running Simulations

1. Select a molecule from your saved collection
2. Click the simulation icon to configure simulation parameters
3. Choose between docking, molecular dynamics, or ADMET prediction
4. Review simulation results and make data-driven refinements

### Collaboration

1. Create a new project from the Projects section
2. Add team members to collaborate
3. Share molecules, comments, and simulation results
4. Track all changes through the activity feed

## API Documentation

The backend provides a comprehensive REST API:

- `/api/auth`: Authentication and user management
- `/api/ai`: Claude AI integration for molecule generation 
- `/api/drug-design`: Molecular design and analysis
- `/api/simulation`: Computational chemistry simulations
- `/api/chemical-data`: ChEMBL database integration
- `/api/literature`: PubMed search and analysis
- `/api/regulatory`: Regulatory pathway analysis
- `/api/collaboration`: Project management and team collaboration

For detailed API documentation, refer to the API guide in the docs folder.

## Technologies Used

- **Frontend**: React, Material-UI, 3Dmol.js
- **Backend**: Node.js, Express
- **AI**: Claude by Anthropic
- **Cheminformatics**: RDKit
- **Authentication**: JWT
- **Database**: File-based JSON storage (production version would use MongoDB/PostgreSQL)
- **External APIs**: PubMed, ChEMBL

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Anthropic for providing the Claude API
- RDKit project for the cheminformatics toolkit
- 3Dmol.js for molecular visualization
- ChEMBL and PubMed for providing chemical and literature databases 