const express = require('express');
const router = express.Router();
const fs = require('fs');
const path = require('path');
const { v4: uuidv4 } = require('uuid');
const axios = require('axios');
const logger = require('../utils/logger');

// Create molecules directory if it doesn't exist
const moleculesDir = path.join(__dirname, '../data/molecules');
if (!fs.existsSync(moleculesDir)) {
  fs.mkdirSync(moleculesDir, { recursive: true });
}

// Helper functions for molecule management
const getMoleculeFilePath = (id) => path.join(moleculesDir, `${id}.json`);

const saveMolecule = (molecule) => {
  // Generate ID if not provided
  const id = molecule.id || uuidv4();
  const newMolecule = {
    ...molecule,
    id,
    dateCreated: molecule.dateCreated || new Date().toISOString()
  };
  
  // Save to file
  const moleculeFile = getMoleculeFilePath(id);
  fs.writeFileSync(moleculeFile, JSON.stringify(newMolecule, null, 2));
  
  return newMolecule;
};

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

// Generate molecule based on requirements
router.post('/generate', async (req, res) => {
  try {
    const { requirements, targetReceptor } = req.body;
    
    if (!requirements) {
      return res.status(400).json({ error: 'Design requirements are required' });
    }
    
    // Call the Claude API through our AI service to generate molecules
    const claudeResponse = await axios.post('http://localhost:5000/api/ai/generate-molecule', {
      requirements
    });
    
    // Get the generated molecules
    const generatedMolecules = claudeResponse.data.molecules || [];
    
    // If target receptor is specified, run docking for each molecule
    if (targetReceptor && generatedMolecules.length > 0) {
      // Call the AI service to analyze and dock the molecules
      const analysisResponse = await axios.post('http://localhost:5000/api/ai/analyze-molecules', {
        molecules: generatedMolecules,
        targetReceptor
      });
      
      res.json({
        requestId: claudeResponse.data.requestId,
        molecules: generatedMolecules,
        analysis: analysisResponse.data,
        thinkingProcess: claudeResponse.data.thinking
      });
    } else {
      // Return just the generated molecules
      res.json({
        requestId: claudeResponse.data.requestId,
        molecules: generatedMolecules,
        thinkingProcess: claudeResponse.data.thinking
      });
    }
  } catch (error) {
    logger.error(`Error generating molecule: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Refine an existing molecule
router.post('/refine', async (req, res) => {
  try {
    const { smiles, requirements, modifications } = req.body;
    
    if (!smiles || !requirements) {
      return res.status(400).json({ 
        error: 'SMILES structure and refinement requirements are required' 
      });
    }
    
    // Build the refinement prompt
    let refinementPrompt = `Refine the following molecule (SMILES: ${smiles}) according to these requirements: ${requirements}`;
    
    if (modifications) {
      refinementPrompt += `\n\nSpecific modifications to make: ${modifications}`;
    }
    
    // Call Claude API through our AI service
    const claudeResponse = await axios.post('http://localhost:5000/api/ai/ask', {
      question: refinementPrompt,
      context: `Original SMILES: ${smiles}`
    });
    
    // Extract SMILES from Claude's response
    const responseText = claudeResponse.data.response;
    
    // Basic regex to extract SMILES from text
    const smilesRegex = /\b([A-Za-z0-9@+\-\[\]\(\)\\\/%=#$!.~{},*]+)\b/g;
    const matches = responseText.match(smilesRegex) || [];
    
    // Filter likely SMILES strings (basic heuristic)
    const potentialSmiles = matches.filter(match => 
      match.length > 5 && 
      !match.includes('http') && 
      (match.includes('C') || match.includes('N')) && 
      (match.includes('(') || match.includes('='))
    );
    
    // Validate the potential SMILES structures
    const validatedMolecules = [];
    
    for (const potentialSmile of potentialSmiles) {
      try {
        // Validate with RDKit via simulation API
        const validationResponse = await axios.post('http://localhost:5000/api/simulation/properties', {
          smiles: potentialSmile
        });
        
        if (!validationResponse.data.error) {
          // Valid SMILES
          validatedMolecules.push({
            smiles: potentialSmile,
            properties: validationResponse.data
          });
        }
      } catch (err) {
        // Invalid SMILES, skip
        logger.error(`Invalid SMILES: ${potentialSmile}`);
      }
    }
    
    res.json({
      originalSmiles: smiles,
      refinedMolecules: validatedMolecules,
      claudeResponse: responseText
    });
  } catch (error) {
    logger.error(`Error refining molecule: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Regulatory analysis for a molecule
router.post('/regulatory-analysis', async (req, res) => {
  try {
    const { smiles, targetIndication } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES structure is required' });
    }
    
    // Call ADMET prediction API
    const admetResponse = await axios.post('http://localhost:5000/api/simulation/admet', {
      smiles
    });
    
    // Get properties
    const propertiesResponse = await axios.post('http://localhost:5000/api/simulation/properties', {
      smiles
    });
    
    // Ask Claude for regulatory analysis
    const systemPrompt = `You are an expert in pharmaceutical regulatory affairs and drug development. Provide a detailed analysis of the regulatory pathway for a candidate molecule based on its properties and ADMET profile.`;
    
    const userPrompt = `Analyze the regulatory pathway for this molecule:
SMILES: ${smiles}
Properties: ${JSON.stringify(propertiesResponse.data)}
ADMET Profile: ${JSON.stringify(admetResponse.data)}
${targetIndication ? `Target Indication: ${targetIndication}` : ''}

Provide:
1. Potential regulatory challenges based on the molecule's properties
2. Required studies for IND submission
3. Estimated timeline for clinical development
4. Potential for expedited pathways
5. Key concerns that might arise during regulatory review`;
    
    // Call Claude through our AI service
    const claudeResponse = await axios.post('http://localhost:5000/api/ai/ask', {
      question: userPrompt
    });
    
    res.json({
      smiles,
      properties: propertiesResponse.data,
      admet: admetResponse.data,
      regulatoryAnalysis: claudeResponse.data.response
    });
  } catch (error) {
    logger.error(`Error performing regulatory analysis: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// CRUD operations for molecules
// Get all saved molecules
router.get('/molecules', (req, res) => {
  try {
    const molecules = getAllMolecules();
    res.json(molecules);
  } catch (err) {
    logger.error(`Error getting all molecules: ${err.message}`);
    res.status(500).json({ error: 'Failed to retrieve molecules' });
  }
});

// Get specific molecule by ID
router.get('/molecules/:id', (req, res) => {
  try {
    const { id } = req.params;
    const moleculeFile = getMoleculeFilePath(id);
    
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

// Save a molecule
router.post('/molecules', async (req, res) => {
  try {
    const moleculeData = req.body;
    
    if (!moleculeData.smiles) {
      return res.status(400).json({ error: 'SMILES structure is required' });
    }
    
    // If no properties are provided, calculate them
    if (!moleculeData.properties) {
      try {
        const propertiesResponse = await axios.post('http://localhost:5000/api/simulation/properties', {
          smiles: moleculeData.smiles
        });
        
        if (!propertiesResponse.data.error) {
          moleculeData.properties = propertiesResponse.data;
        }
      } catch (err) {
        logger.error(`Error calculating properties: ${err.message}`);
        // Continue without properties if calculation fails
      }
    }
    
    const savedMolecule = saveMolecule(moleculeData);
    res.status(201).json(savedMolecule);
  } catch (err) {
    logger.error(`Error saving molecule: ${err.message}`);
    res.status(500).json({ error: 'Failed to save molecule' });
  }
});

// Update molecule
router.put('/molecules/:id', (req, res) => {
  try {
    const { id } = req.params;
    const moleculeData = req.body;
    const moleculeFile = getMoleculeFilePath(id);
    
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
router.delete('/molecules/:id', (req, res) => {
  try {
    const { id } = req.params;
    const moleculeFile = getMoleculeFilePath(id);
    
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

module.exports = router; 