const express = require('express');
const cors = require('cors');
const morgan = require('morgan');
const helmet = require('helmet');
const dotenv = require('dotenv');
const { spawn } = require('child_process');

// Load environment variables
dotenv.config();

// Set up Express app
const app = express();
const PORT = process.env.PORT || 5000;

// Middleware
app.use(cors());
app.use(helmet());
app.use(morgan('dev'));
app.use(express.json());

// Import API routes
const literatureRoutes = require('./api/literature');
const chemicalDataRoutes = require('./api/chemicalData');
const drugDesignRoutes = require('./api/drugDesign');
const simulationRoutes = require('./api/simulation');
const regulatoryRoutes = require('./api/regulatory');
const aiRoutes = require('./api/ai');
const { router: authRoutes } = require('./api/auth');
const collaborationRoutes = require('./api/collaboration');
const moleculeRoutes = require('./routes/molecule');
const similarityRoutes = require('./api/similarity');

// API Routes
app.use('/api/literature', literatureRoutes);
app.use('/api/chemical-data', chemicalDataRoutes);
app.use('/api/drug-design', drugDesignRoutes);
app.use('/api/simulation', simulationRoutes);
app.use('/api/regulatory', regulatoryRoutes);
app.use('/api/ai', aiRoutes);
app.use('/api/auth', authRoutes);
app.use('/api/collaboration', collaborationRoutes);
app.use('/api/molecule', moleculeRoutes);
app.use('/api/similarity', similarityRoutes);

// Health check endpoint
app.get('/api/health', (req, res) => {
  res.status(200).json({ status: 'up', message: 'Breaking Good API is running' });
});

// RDKit Python bridge function
const runRDKitScript = (scriptName, args) => {
  return new Promise((resolve, reject) => {
    const pythonProcess = spawn('python', [`./utils/rdkit/${scriptName}.py`, ...args]);
    
    let result = '';
    let error = '';
    
    pythonProcess.stdout.on('data', (data) => {
      result += data.toString();
    });
    
    pythonProcess.stderr.on('data', (data) => {
      error += data.toString();
    });
    
    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        reject(new Error(`RDKit script exited with code ${code}: ${error}`));
      } else {
        try {
          resolve(JSON.parse(result));
        } catch (e) {
          resolve(result);
        }
      }
    });
  });
};

// Make RDKit bridge available to other modules
app.locals.runRDKitScript = runRDKitScript;

// Start server with better error handling
const server = app.listen(PORT, () => {
  console.log(`Breaking Good API running on port ${PORT}`);
});

server.on('error', (e) => {
  if (e.code === 'EADDRINUSE') {
    console.error(`Port ${PORT} is already in use. Please close the other application using this port or change the PORT in .env file.`);
    process.exit(1);
  } else {
    console.error('Server error:', e);
  }
});

module.exports = app; 