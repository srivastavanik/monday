const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const os = require('os');
const { v4: uuidv4 } = require('uuid');

// Middleware to validate request body
const validateSimilarityRequest = (req, res, next) => {
  const { query, targets, fingerprint, metric, threshold } = req.body;
  
  // Check if query molecule is provided
  if (!query) {
    return res.status(400).json({ error: 'Query molecule SMILES string is required' });
  }
  
  // Check if targets are provided correctly
  if (!targets || !Array.isArray(targets)) {
    return res.status(400).json({ error: 'Target molecules must be provided as an array' });
  }
  
  // Validate optional parameters
  if (fingerprint && !['morgan', 'maccs', 'rdkit', 'avalon'].includes(fingerprint)) {
    return res.status(400).json({ error: 'Invalid fingerprint type' });
  }
  
  if (metric && !['tanimoto', 'dice', 'cosine', 'sokal', 'russel'].includes(metric)) {
    return res.status(400).json({ error: 'Invalid similarity metric' });
  }
  
  if (threshold && (threshold < 0 || threshold > 1)) {
    return res.status(400).json({ error: 'Threshold must be between 0 and 1' });
  }
  
  next();
};

// API endpoint for similarity search
router.post('/search', validateSimilarityRequest, async (req, res) => {
  try {
    const { query, targets, fingerprint = 'morgan', metric = 'tanimoto', threshold = 0.7, maxResults = 50 } = req.body;
    
    // Create a temporary JSON file to store target compounds
    const tempDir = os.tmpdir();
    const tempFilePath = path.join(tempDir, `targets-${uuidv4()}.json`);
    
    // Format targets as expected by the Python script
    const formattedTargets = targets.map(target => {
      // Ensure each target has at least a SMILES string
      if (typeof target === 'string') {
        return { smiles: target };
      } else if (typeof target === 'object' && target.smiles) {
        return target;
      } else {
        return null;
      }
    }).filter(Boolean);
    
    // Write targets to temp file
    fs.writeFileSync(tempFilePath, JSON.stringify(formattedTargets));
    
    // Build arguments for the Python script
    const scriptPath = path.join(__dirname, '../utils/rdkit/similarity_search.py');
    const args = [
      scriptPath,
      '--query', query,
      '--target_file', tempFilePath,
      '--fingerprint', fingerprint,
      '--metric', metric,
      '--threshold', threshold.toString(),
      '--max_results', maxResults.toString()
    ];
    
    // Execute the Python script
    const pythonProcess = spawn('python', args);
    
    let resultData = '';
    let errorData = '';
    
    pythonProcess.stdout.on('data', (data) => {
      resultData += data.toString();
    });
    
    pythonProcess.stderr.on('data', (data) => {
      errorData += data.toString();
    });
    
    pythonProcess.on('close', (code) => {
      // Clean up temporary file
      try {
        fs.unlinkSync(tempFilePath);
      } catch (err) {
        console.error('Error deleting temporary file:', err);
      }
      
      if (code !== 0) {
        return res.status(500).json({
          error: 'Similarity search failed',
          details: errorData
        });
      }
      
      try {
        const results = JSON.parse(resultData);
        res.json(results);
      } catch (err) {
        res.status(500).json({
          error: 'Failed to parse similarity search results',
          details: errorData || err.message
        });
      }
    });
  } catch (err) {
    res.status(500).json({
      error: 'Server error',
      details: err.message
    });
  }
});

// API endpoint for multi-reference similarity search
router.post('/multi-reference', validateSimilarityRequest, async (req, res) => {
  try {
    const { 
      queries, 
      targets, 
      fingerprint = 'morgan', 
      metric = 'tanimoto', 
      threshold = 0.7, 
      maxResults = 50,
      aggregation = 'max'
    } = req.body;
    
    // Create a temporary file to store the Python script input
    const tempDir = os.tmpdir();
    const inputFilePath = path.join(tempDir, `input-${uuidv4()}.json`);
    const outputFilePath = path.join(tempDir, `output-${uuidv4()}.json`);
    
    // Format targets as expected by the Python script
    const formattedTargets = targets.map(target => {
      if (typeof target === 'string') {
        return { smiles: target };
      } else if (typeof target === 'object' && target.smiles) {
        return target;
      } else {
        return null;
      }
    }).filter(Boolean);
    
    // Write input data for the Python script
    const inputData = {
      queries: Array.isArray(queries) ? queries : [queries],
      targets: formattedTargets,
      fingerprint,
      metric,
      threshold,
      maxResults,
      aggregation
    };
    
    fs.writeFileSync(inputFilePath, JSON.stringify(inputData));
    
    // Create a Python script to use the MolecularSimilaritySearch class
    const pythonScript = `
import json
import sys
sys.path.append('${path.join(__dirname, '../utils/rdkit')}')
from similarity_search import MolecularSimilaritySearch

# Load input data
with open('${inputFilePath}', 'r') as f:
    input_data = json.load(f)

# Initialize the similarity search
sim_search = MolecularSimilaritySearch(
    fingerprint_type=input_data['fingerprint'],
    similarity_metric=input_data['metric']
)

# Perform multi-reference search
results = sim_search.search_multiple_references(
    input_data['queries'],
    input_data['targets'],
    threshold=input_data['threshold'],
    max_results=input_data['maxResults'],
    aggregation=input_data['aggregation']
)

# Write results
with open('${outputFilePath}', 'w') as f:
    json.dump(results, f)
    `;
    
    const scriptPath = path.join(tempDir, `script-${uuidv4()}.py`);
    fs.writeFileSync(scriptPath, pythonScript);
    
    // Execute the Python script
    const pythonProcess = spawn('python', [scriptPath]);
    
    let errorData = '';
    
    pythonProcess.stderr.on('data', (data) => {
      errorData += data.toString();
    });
    
    pythonProcess.on('close', (code) => {
      // Clean up temporary files
      try {
        fs.unlinkSync(inputFilePath);
        fs.unlinkSync(scriptPath);
      } catch (err) {
        console.error('Error deleting temporary files:', err);
      }
      
      if (code !== 0) {
        return res.status(500).json({
          error: 'Multi-reference similarity search failed',
          details: errorData
        });
      }
      
      try {
        // Read the results
        const resultData = fs.readFileSync(outputFilePath, 'utf8');
        fs.unlinkSync(outputFilePath);
        
        const results = JSON.parse(resultData);
        res.json(results);
      } catch (err) {
        res.status(500).json({
          error: 'Failed to parse similarity search results',
          details: errorData || err.message
        });
      }
    });
  } catch (err) {
    res.status(500).json({
      error: 'Server error',
      details: err.message
    });
  }
});

// API endpoint for diversity selection
router.post('/diversity', async (req, res) => {
  try {
    const { 
      molecules, 
      numPicks = 10, 
      fingerprint = 'morgan', 
      metric = 'tanimoto' 
    } = req.body;
    
    if (!molecules || !Array.isArray(molecules) || molecules.length === 0) {
      return res.status(400).json({ error: 'A list of molecule SMILES strings is required' });
    }
    
    // Create a temporary file to store the Python script input
    const tempDir = os.tmpdir();
    const inputFilePath = path.join(tempDir, `input-${uuidv4()}.json`);
    const outputFilePath = path.join(tempDir, `output-${uuidv4()}.json`);
    
    // Format molecules
    const formattedMolecules = molecules.map(mol => {
      if (typeof mol === 'string') {
        return mol;
      } else if (typeof mol === 'object' && mol.smiles) {
        return mol.smiles;
      } else {
        return null;
      }
    }).filter(Boolean);
    
    // Write input data for the Python script
    const inputData = {
      molecules: formattedMolecules,
      numPicks: numPicks,
      fingerprint,
      metric
    };
    
    fs.writeFileSync(inputFilePath, JSON.stringify(inputData));
    
    // Create a Python script to use the MolecularSimilaritySearch class
    const pythonScript = `
import json
import sys
sys.path.append('${path.join(__dirname, '../utils/rdkit')}')
from similarity_search import MolecularSimilaritySearch

# Load input data
with open('${inputFilePath}', 'r') as f:
    input_data = json.load(f)

# Initialize the similarity search
sim_search = MolecularSimilaritySearch(
    fingerprint_type=input_data['fingerprint'],
    similarity_metric=input_data['metric']
)

# Perform diversity selection
indices = sim_search.diversity_selection(
    input_data['molecules'],
    num_picks=input_data['numPicks']
)

# Create results with selected molecules
if isinstance(indices, dict) and 'error' in indices:
    results = indices
else:
    results = {
        'selected_indices': indices,
        'selected_molecules': [input_data['molecules'][i] for i in indices]
    }

# Write results
with open('${outputFilePath}', 'w') as f:
    json.dump(results, f)
    `;
    
    const scriptPath = path.join(tempDir, `script-${uuidv4()}.py`);
    fs.writeFileSync(scriptPath, pythonScript);
    
    // Execute the Python script
    const pythonProcess = spawn('python', [scriptPath]);
    
    let errorData = '';
    
    pythonProcess.stderr.on('data', (data) => {
      errorData += data.toString();
    });
    
    pythonProcess.on('close', (code) => {
      // Clean up temporary files
      try {
        fs.unlinkSync(inputFilePath);
        fs.unlinkSync(scriptPath);
      } catch (err) {
        console.error('Error deleting temporary files:', err);
      }
      
      if (code !== 0) {
        return res.status(500).json({
          error: 'Diversity selection failed',
          details: errorData
        });
      }
      
      try {
        // Read the results
        const resultData = fs.readFileSync(outputFilePath, 'utf8');
        fs.unlinkSync(outputFilePath);
        
        const results = JSON.parse(resultData);
        res.json(results);
      } catch (err) {
        res.status(500).json({
          error: 'Failed to parse diversity selection results',
          details: errorData || err.message
        });
      }
    });
  } catch (err) {
    res.status(500).json({
      error: 'Server error',
      details: err.message
    });
  }
});

module.exports = router; 