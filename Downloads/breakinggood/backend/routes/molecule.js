const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { v4: uuidv4 } = require('uuid');
const logger = require('../utils/logger');

/**
 * @route   POST /api/molecule/descriptors
 * @desc    Calculate molecular descriptors for a given SMILES string
 * @access  Public
 */
router.post('/descriptors', async (req, res) => {
  try {
    const { smiles } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES string is required' });
    }
    
    // Validate SMILES string (basic validation)
    if (!/^[A-Za-z0-9@+\-\[\]\(\)\\\/%=#$:.~!,]+$/.test(smiles)) {
      return res.status(400).json({ error: 'Invalid SMILES string format' });
    }
    
    // Path to the Python script
    const scriptPath = path.join(__dirname, '../utils/rdkit/descriptors.py');
    
    // Ensure the script exists
    if (!fs.existsSync(scriptPath)) {
      return res.status(500).json({ error: 'RDKit descriptor calculator not found' });
    }
    
    // Spawn a Python process to calculate descriptors
    const pythonProcess = spawn('python', [scriptPath, smiles]);
    
    let result = '';
    let errorOutput = '';
    
    // Collect data from stdout
    pythonProcess.stdout.on('data', (data) => {
      result += data.toString();
    });
    
    // Collect error data from stderr
    pythonProcess.stderr.on('data', (data) => {
      errorOutput += data.toString();
    });
    
    // Handle process completion
    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        console.error(`Python process exited with code ${code}`);
        console.error(`Error output: ${errorOutput}`);
        return res.status(500).json({ 
          error: 'Error calculating molecular descriptors',
          details: errorOutput
        });
      }
      
      try {
        // Parse the JSON result
        const descriptors = JSON.parse(result);
        
        if (descriptors.error) {
          return res.status(400).json({ error: descriptors.error });
        }
        
        return res.json(descriptors);
      } catch (parseError) {
        console.error('Error parsing Python output:', parseError);
        return res.status(500).json({ 
          error: 'Error parsing descriptor results',
          details: result
        });
      }
    });
    
  } catch (error) {
    console.error('Server error:', error);
    res.status(500).json({ error: 'Server error' });
  }
});

/**
 * @route   POST /api/molecule/generate
 * @desc    Generate molecules based on input parameters
 * @access  Public
 */
router.post('/generate', async (req, res) => {
  try {
    const { 
      baseStructure,
      targetReceptor,
      similarityThreshold,
      maxMolecules,
      optimizationTarget
    } = req.body;
    
    // In a real implementation, this would call a molecule generation service
    // For now, we'll just return some dummy molecules as a placeholder
    
    // Simulate a delay for processing
    await new Promise(resolve => setTimeout(resolve, 2000));
    
    // Example response with generated molecules
    const generatedMolecules = [
      {
        id: 'mol1',
        smiles: 'CC1=C(C(=O)N(C)C2=CC=CC=C21)SC3=CC=CC=C3',
        name: 'Generated Compound 1',
        similarity: 0.85,
        predictedActivity: 9.2,
        properties: {
          molecularWeight: 308.4,
          logP: 3.8,
          hBondDonors: 0,
          hBondAcceptors: 2,
          rotableBonds: 2,
          qed: 0.72
        }
      },
      {
        id: 'mol2',
        smiles: 'CC1=NN(C(=O)C1)C2=CC=C(C=C2)S(=O)(=O)N',
        name: 'Generated Compound 2',
        similarity: 0.78,
        predictedActivity: 8.7,
        properties: {
          molecularWeight: 279.3,
          logP: 1.2,
          hBondDonors: 1,
          hBondAcceptors: 6,
          rotableBonds: 3,
          qed: 0.65
        }
      },
      {
        id: 'mol3',
        smiles: 'C1=CC=C(C=C1)CN2C=C(N=N2)C(=O)NC3=CC=CC=C3',
        name: 'Generated Compound 3',
        similarity: 0.92,
        predictedActivity: 9.8,
        properties: {
          molecularWeight: 290.3,
          logP: 2.9,
          hBondDonors: 1,
          hBondAcceptors: 4,
          rotableBonds: 5,
          qed: 0.81
        }
      }
    ];
    
    res.json({ 
      molecules: generatedMolecules,
      message: 'Successfully generated molecules'
    });
    
  } catch (error) {
    console.error('Server error:', error);
    res.status(500).json({ error: 'Server error' });
  }
});

// Create molecules directory if it doesn't exist
const moleculesDir = path.join(__dirname, '../data/molecules');
if (!fs.existsSync(moleculesDir)) {
  fs.mkdirSync(moleculesDir, { recursive: true });
}

// Helper to get all molecules
const getAllMolecules = () => {
  try {
    const files = fs.readdirSync(moleculesDir);
    const molecules = files
      .filter(file => file.endsWith('.json'))
      .map(file => {
        try {
          const data = fs.readFileSync(path.join(moleculesDir, file), 'utf8');
          return JSON.parse(data);
        } catch (err) {
          logger.error(`Error reading molecule file ${file}: ${err.message}`);
          return null;
        }
      })
      .filter(molecule => molecule !== null);
    
    return molecules;
  } catch (err) {
    logger.error(`Error reading molecules directory: ${err.message}`);
    return [];
  }
};

// Get all molecules
router.get('/', (req, res) => {
  try {
    const molecules = getAllMolecules();
    res.json(molecules);
  } catch (err) {
    logger.error(`Error getting all molecules: ${err.message}`);
    res.status(500).json({ error: 'Failed to retrieve molecules' });
  }
});

// Get molecule by ID
router.get('/:id', (req, res) => {
  try {
    const { id } = req.params;
    const moleculeFile = path.join(moleculesDir, `${id}.json`);
    
    if (!fs.existsSync(moleculeFile)) {
      return res.status(404).json({ error: 'Molecule not found' });
    }
    
    const moleculeData = JSON.parse(fs.readFileSync(moleculeFile, 'utf8'));
    res.json(moleculeData);
  } catch (err) {
    logger.error(`Error getting molecule: ${err.message}`);
    res.status(500).json({ error: 'Failed to retrieve molecule' });
  }
});

// Create new molecule
router.post('/', (req, res) => {
  try {
    const moleculeData = req.body;
    
    if (!moleculeData.smiles) {
      return res.status(400).json({ error: 'SMILES structure is required' });
    }
    
    // Generate ID if not provided
    const id = moleculeData.id || uuidv4();
    const newMolecule = {
      ...moleculeData,
      id,
      dateCreated: moleculeData.dateCreated || new Date().toISOString()
    };
    
    // Save to file
    const moleculeFile = path.join(moleculesDir, `${id}.json`);
    fs.writeFileSync(moleculeFile, JSON.stringify(newMolecule, null, 2));
    
    res.status(201).json(newMolecule);
  } catch (err) {
    logger.error(`Error creating molecule: ${err.message}`);
    res.status(500).json({ error: 'Failed to create molecule' });
  }
});

// Update molecule
router.put('/:id', (req, res) => {
  try {
    const { id } = req.params;
    const moleculeData = req.body;
    const moleculeFile = path.join(moleculesDir, `${id}.json`);
    
    if (!fs.existsSync(moleculeFile)) {
      return res.status(404).json({ error: 'Molecule not found' });
    }
    
    const existingMolecule = JSON.parse(fs.readFileSync(moleculeFile, 'utf8'));
    
    const updatedMolecule = {
      ...existingMolecule,
      ...moleculeData,
      id, // Ensure ID remains the same
      dateModified: new Date().toISOString()
    };
    
    fs.writeFileSync(moleculeFile, JSON.stringify(updatedMolecule, null, 2));
    
    res.json(updatedMolecule);
  } catch (err) {
    logger.error(`Error updating molecule: ${err.message}`);
    res.status(500).json({ error: 'Failed to update molecule' });
  }
});

// Delete molecule
router.delete('/:id', (req, res) => {
  try {
    const { id } = req.params;
    const moleculeFile = path.join(moleculesDir, `${id}.json`);
    
    if (!fs.existsSync(moleculeFile)) {
      return res.status(404).json({ error: 'Molecule not found' });
    }
    
    fs.unlinkSync(moleculeFile);
    
    res.json({ message: 'Molecule deleted successfully' });
  } catch (err) {
    logger.error(`Error deleting molecule: ${err.message}`);
    res.status(500).json({ error: 'Failed to delete molecule' });
  }
});

// Search molecules
router.get('/search/:query', (req, res) => {
  try {
    const { query } = req.params;
    const molecules = getAllMolecules();
    
    const matchingMolecules = molecules.filter(molecule => {
      return (
        (molecule.name && molecule.name.toLowerCase().includes(query.toLowerCase())) ||
        (molecule.smiles && molecule.smiles.toLowerCase().includes(query.toLowerCase())) ||
        (molecule.description && molecule.description.toLowerCase().includes(query.toLowerCase()))
      );
    });
    
    res.json(matchingMolecules);
  } catch (err) {
    logger.error(`Error searching molecules: ${err.message}`);
    res.status(500).json({ error: 'Failed to search molecules' });
  }
});

// Batch import molecules
router.post('/batch', (req, res) => {
  try {
    const { molecules } = req.body;
    
    if (!Array.isArray(molecules) || molecules.length === 0) {
      return res.status(400).json({ error: 'Valid molecules array is required' });
    }
    
    const savedMolecules = [];
    
    for (const moleculeData of molecules) {
      if (!moleculeData.smiles) {
        continue; // Skip molecules without SMILES
      }
      
      // Generate ID if not provided
      const id = moleculeData.id || uuidv4();
      const newMolecule = {
        ...moleculeData,
        id,
        dateCreated: moleculeData.dateCreated || new Date().toISOString()
      };
      
      // Save to file
      const moleculeFile = path.join(moleculesDir, `${id}.json`);
      fs.writeFileSync(moleculeFile, JSON.stringify(newMolecule, null, 2));
      
      savedMolecules.push(newMolecule);
    }
    
    res.status(201).json({
      message: `Successfully imported ${savedMolecules.length} molecules`,
      molecules: savedMolecules
    });
  } catch (err) {
    logger.error(`Error batch importing molecules: ${err.message}`);
    res.status(500).json({ error: 'Failed to import molecules' });
  }
});

module.exports = router; 