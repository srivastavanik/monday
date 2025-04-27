const express = require('express');
const router = express.Router();
const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const { v4: uuidv4 } = require('uuid');
const logger = require('../utils/logger');

// Create results directory if it doesn't exist
const resultsDir = path.join(__dirname, '../data/simulation_results');
if (!fs.existsSync(resultsDir)) {
  fs.mkdirSync(resultsDir, { recursive: true });
}

// Path to RDKit scripts
const MOLECULAR_PROPERTIES_SCRIPT = path.join(__dirname, '../utils/rdkit/molecular_properties.py');
const ADMET_PREDICTION_SCRIPT = path.join(__dirname, '../utils/rdkit/admet_prediction.py');
const MOLECULAR_DOCKING_SCRIPT = path.join(__dirname, '../utils/rdkit/molecular_docking.py');
const SIMILARITY_SEARCH_SCRIPT = path.join(__dirname, '../utils/rdkit/similarity_search.py');
const MOLECULE_OPERATIONS_SCRIPT = path.join(__dirname, '../utils/rdkit/molecule_operations.py');
const STRUCTURE_GENERATION_SCRIPT = path.join(__dirname, '../utils/rdkit/structure_generation.py');
const FORMAT_CONVERTER_SCRIPT = path.join(__dirname, '../utils/rdkit/format_converter.py');

// Helper function to create a temporary file
const createTempFile = (prefix, suffix) => {
  const tempId = uuidv4();
  return path.join(resultsDir, `${prefix}_${tempId}${suffix}`);
};

// Helper function to run RDKit Python scripts
const runRDKitScript = (scriptPath, args = []) => {
  return new Promise((resolve, reject) => {
    const pythonProcess = spawn('python', [scriptPath, ...args]);
    
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

// Calculate molecular properties
router.post('/properties', async (req, res) => {
  try {
    const { smiles } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES string is required' });
    }
    
    // Run molecular properties calculation script
    const properties = await runRDKitScript(MOLECULAR_PROPERTIES_SCRIPT, [smiles]);
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `properties_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({ 
      id: resultId,
      timestamp: new Date().toISOString(),
      smiles,
      properties
    }, null, 2));
    
    return res.json({
      id: resultId,
      properties
    });
    
  } catch (error) {
    logger.error(`Error calculating molecular properties: ${error.message}`);
    return res.status(500).json({
      error: 'Error calculating molecular properties',
      details: error.message
    });
  }
});

// Predict ADMET properties
router.post('/admet', async (req, res) => {
  try {
    const { smiles } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES string is required' });
    }
    
    // Run ADMET prediction script
    const admetResults = await runRDKitScript(ADMET_PREDICTION_SCRIPT, [smiles]);
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `admet_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({
      id: resultId,
      timestamp: new Date().toISOString(),
      smiles,
      admet: admetResults
    }, null, 2));
    
    return res.json({
      id: resultId,
      admet: admetResults
    });
    
  } catch (error) {
    logger.error(`Error predicting ADMET properties: ${error.message}`);
    return res.status(500).json({
      error: 'Error predicting ADMET properties',
      details: error.message
    });
  }
});

// Run molecular docking
router.post('/docking', async (req, res) => {
  try {
    const { 
      receptorPdb, 
      ligandSmiles, 
      exhaustiveness = 8, 
      centerX = 0, 
      centerY = 0, 
      centerZ = 0, 
      sizeX = 20, 
      sizeY = 20, 
      sizeZ = 20 
    } = req.body;
    
    if (!receptorPdb) {
      return res.status(400).json({ error: 'Receptor PDB data is required' });
    }
    
    if (!ligandSmiles) {
      return res.status(400).json({ error: 'Ligand SMILES string is required' });
    }
    
    // Save receptor PDB to temporary file if it's not a path
    const tempId = uuidv4();
    let receptorPath = receptorPdb;
    if (!fs.existsSync(receptorPdb)) {
      const receptorFile = path.join(resultsDir, `receptor_${tempId}.pdb`);
      fs.writeFileSync(receptorFile, receptorPdb);
      receptorPath = receptorFile;
    }
    
    // Run molecular docking script
    const dockingResults = await runRDKitScript(
      MOLECULAR_DOCKING_SCRIPT, 
      [
        ligandSmiles,
        receptorPath,
        exhaustiveness.toString(),
        centerX.toString(),
        centerY.toString(),
        centerZ.toString(),
        sizeX.toString(),
        sizeY.toString(),
        sizeZ.toString()
      ]
    );
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `docking_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({
      id: resultId,
      timestamp: new Date().toISOString(),
      receptorPath,
      ligandSmiles,
      params: {
        exhaustiveness,
        center: [centerX, centerY, centerZ],
        size: [sizeX, sizeY, sizeZ]
      },
      dockingResults
    }, null, 2));
    
    return res.json({
      id: resultId,
      ...dockingResults
    });
    
  } catch (error) {
    logger.error(`Error running molecular docking: ${error.message}`);
    return res.status(500).json({
      error: 'Error running molecular docking',
      details: error.message
    });
  }
});

// Run similarity search
router.post('/similarity', async (req, res) => {
  try {
    const { querySmiles, threshold = 0.7, limit = 10, database } = req.body;
    
    if (!querySmiles) {
      return res.status(400).json({ error: 'Query SMILES string is required' });
    }
    
    if (!database) {
      return res.status(400).json({ error: 'Database path or molecules array is required' });
    }
    
    // If database is an array of molecules, save to temporary file
    let databasePath = database;
    if (Array.isArray(database)) {
      const dbFile = path.join(resultsDir, `molecules_${uuidv4()}.json`);
      fs.writeFileSync(dbFile, JSON.stringify(database));
      databasePath = dbFile;
    }
    
    // Run similarity search script
    const similarityResults = await runRDKitScript(
      SIMILARITY_SEARCH_SCRIPT, 
      [
        querySmiles, 
        databasePath,
        threshold.toString(),
        limit.toString()
      ]
    );
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `similarity_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({
      id: resultId,
      timestamp: new Date().toISOString(),
      querySmiles,
      threshold,
      limit,
      results: similarityResults
    }, null, 2));
    
    return res.json({
      id: resultId,
      results: similarityResults
    });
    
  } catch (error) {
    logger.error(`Error running similarity search: ${error.message}`);
    return res.status(500).json({
      error: 'Error running similarity search',
      details: error.message
    });
  }
});

// Generate 3D structure
router.post('/generate-3d', async (req, res) => {
  try {
    const { smiles } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES string is required' });
    }
    
    // Run 3D structure generation script
    const structureResults = await runRDKitScript(
      STRUCTURE_GENERATION_SCRIPT, 
      [smiles]
    );
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `structure_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({
      id: resultId,
      timestamp: new Date().toISOString(),
      smiles,
      structureResults
    }, null, 2));
    
    return res.json({
      id: resultId,
      structureResults
    });
    
  } catch (error) {
    logger.error(`Error generating 3D structure: ${error.message}`);
    return res.status(500).json({
      error: 'Error generating 3D structure',
      details: error.message
    });
  }
});

// Convert between molecular formats
router.post('/convert', async (req, res) => {
  const { input, inputFormat, outputFormat } = req.body;
  try {
    const molInputPreview = (typeof input === 'string') ? input.substring(0, 200) + (input.length > 200 ? '...' : '') : '[Input not a string]';
    logger.debug(`Received /convert request. Format: ${inputFormat} -> ${outputFormat}. Preview: ${molInputPreview}`);

    if (!input) {
      return res.status(400).json({ error: 'Input molecule data is required' });
    }
    if (!inputFormat || !outputFormat) {
      return res.status(400).json({ error: 'Input and output formats are required' });
    }
    
    // --- Use the new FORMAT_CONVERTER_SCRIPT --- 
    // Script expects positional args: input_data, input_format, output_format
    const scriptArgs = [input, inputFormat, outputFormat];
    logger.info(`Running format converter script...`);
    logger.debug('Script args:', scriptArgs);

    const conversionResult = await runRDKitScript(
      FORMAT_CONVERTER_SCRIPT, // Call the correct script
      scriptArgs
    );
    // -------------------------------------------

    // --- Check for error reported by the script itself --- 
    if (conversionResult && conversionResult.error) {
        logger.error(`Format conversion script reported error: ${conversionResult.error}`);
        // Return 400 for bad input format/data, 500 for others?
        // Let's assume script errors are due to bad input for now.
        return res.status(400).json({ 
            error: 'Format conversion failed', 
            details: conversionResult.error 
        });
    }
    // ----------------------------------------------------
    
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `conversion_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({
      id: resultId,
      timestamp: new Date().toISOString(),
      inputData: input, 
      inputFormat,
      outputFormat,
      // Use the output field from the script's JSON response
      result: conversionResult.output 
    }, null, 2));
    
    return res.json({
      id: resultId,
      // --- Return the 'output' field from the script's JSON --- 
      output: conversionResult.output 
    });
    
  } catch (error) { // Catch errors from runRDKitScript (e.g., script crash)
    const molInputPreview = (typeof input === 'string') ? input.substring(0, 200) + (input.length > 200 ? '...' : '') : '[Input not a string]'; 
    logger.error(`Error running format conversion script (${inputFormat} -> ${outputFormat}):`, {
        errorMessage: error.message,
        inputPreview: molInputPreview,
        stack: error.stack
    });
    const rdkitErrorMatch = error.message.match(/RDKit script exited with code \d+: (.*)/s);
    const specificError = rdkitErrorMatch ? rdkitErrorMatch[1].trim() : 'Internal server error during script execution.';
    
    return res.status(500).json({
      error: 'Format conversion process failed',
      details: specificError
    });
  }
});

// Compare two molecules
router.post('/compare', async (req, res) => {
  try {
    const { smiles1, smiles2, method = 'tanimoto' } = req.body;
    
    if (!smiles1 || !smiles2) {
      return res.status(400).json({ error: 'Two SMILES strings are required for comparison' });
    }
    
    // Run comparison script
    const comparisonResult = await runRDKitScript(
      MOLECULE_OPERATIONS_SCRIPT, 
      ['compare', smiles1, smiles2, method]
    );
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `comparison_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({
      id: resultId,
      timestamp: new Date().toISOString(),
      smiles1,
      smiles2,
      method,
      result: comparisonResult
    }, null, 2));
    
    return res.json({
      id: resultId,
      result: comparisonResult
    });
    
  } catch (error) {
    logger.error(`Error comparing molecules: ${error.message}`);
    return res.status(500).json({
      error: 'Error comparing molecules',
      details: error.message
    });
  }
});

// Run molecular dynamics simulation (Simplified OpenMM version)
router.post('/molecular-dynamics', async (req, res) => {
  try {
    const { 
      smiles, 
      simulationTimePs = 100, // Default to 100 ps
      timeStepFs = 2.0, // Default 2 fs time step
      temperatureK = 300, 
      forceField = 'MMFF94' // Allow specifying force field
    } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES structure is required' });
    }
    
    const simulationId = uuidv4();
    logger.info(`Starting MD simulation for ${smiles} (ID: ${simulationId})`);

    // Calculate steps based on desired time and step size
    const totalSteps = Math.round((simulationTimePs * 1000) / timeStepFs);

    // Call the RDKit script for the short MD simulation
    const result = await runRDKitScript(
        'md', // Use the new 'md' operation
        smiles, 
        totalSteps, // Pass calculated steps
        timeStepFs, 
        temperatureK,
        forceField
    );

    if (result.error) {
      logger.error(`Error in MD simulation (ID: ${simulationId}): ${result.error}`);
      // Include traceback if available
      const details = result.traceback ? `${result.error}\n${result.traceback}` : result.error;
      return res.status(400).json({ error: `MD simulation failed: ${result.error}`, details: details, simulationId });
    }
    
    // Structure the response with results from the simulation
    const simulationData = {
      id: simulationId,
      timestamp: new Date().toISOString(),
      type: 'molecular-dynamics',
      parameters: {
        smiles,
        simulationTimePs,
        timeStepFs,
        temperatureK,
        forceField
      },
      status: 'completed',
      results: { // Nested structure for clarity
        steps_run: result.steps,
        final_potential_energy_kcal_mol: result.final_potential_energy_kcal_mol,
        rmsd_vs_minimized: result.rmsd_vs_minimized,
        // Add other relevant metrics if calculated by the python script
        final_structure_molblock: result.final_molblock // Include the final structure
      }
    };
    
    // Save the results
    const resultsFile = path.join(resultsDir, `${simulationId}.json`);
    fs.writeFileSync(resultsFile, JSON.stringify(simulationData, null, 2));
    logger.info(`MD simulation completed and saved for ID: ${simulationId}`);
    
    // Return the simulation ID and results
    res.json(simulationData);

  } catch (error) {
    logger.error(`Error in molecular dynamics endpoint: ${error.message}`, { stack: error.stack });
    res.status(500).json({ error: `Server error during MD simulation: ${error.message}` });
  }
});

// Get all simulation results
router.get('/results', (req, res) => {
  try {
    const files = fs.readdirSync(resultsDir);
    const results = files.map(file => {
      const filePath = path.join(resultsDir, file);
      const stats = fs.statSync(filePath);
      
      return {
        id: path.basename(file, path.extname(file)),
        file,
        type: path.extname(file).replace('.', ''),
        size: stats.size,
        created: stats.birthtime
      };
    });
    
    return res.json(results);
    
  } catch (error) {
    logger.error(`Error getting simulation results: ${error.message}`);
    return res.status(500).json({
      error: 'Error getting simulation results',
      details: error.message
    });
  }
});

// Get a specific simulation result
router.get('/results/:id', (req, res) => {
  try {
    const { id } = req.params;
    const files = fs.readdirSync(resultsDir);
    
    // Find the file that starts with the ID
    const resultFile = files.find(file => file.startsWith(id) || file.includes(`_${id}`));
    
    if (!resultFile) {
      return res.status(404).json({ error: 'Result not found' });
    }
    
    const filePath = path.join(resultsDir, resultFile);
    const resultData = JSON.parse(fs.readFileSync(filePath, 'utf8'));
    
    return res.json(resultData);
    
  } catch (error) {
    logger.error(`Error getting simulation result: ${error.message}`);
    return res.status(500).json({
      error: 'Error getting simulation result',
      details: error.message
    });
  }
});

// Generate 3D structure from SMILES
router.post('/3d-structure', async (req, res) => {
  try {
    const { smiles } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES string is required' });
    }
    
    // Run structure generation script
    const result = await runRDKitScript(STRUCTURE_GENERATION_SCRIPT, [smiles]);
    
    // Save result
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `3d_structure_${resultId}.json`);
    
    const resultData = {
      id: resultId,
      timestamp: new Date().toISOString(),
      smiles,
      molblock: result.molblock
    };
    
    fs.writeFileSync(resultFile, JSON.stringify(resultData, null, 2));
    
    return res.json({
      id: resultId,
      molblock: result.molblock
    });
    
  } catch (error) {
    logger.error(`Error generating 3D structure: ${error.message}`);
    return res.status(500).json({
      error: 'Error generating 3D structure',
      details: error.message
    });
  }
});

// --- NEW SIMULATE ENDPOINT --- 
router.post('/simulate', async (req, res) => {
  const { smiles, simulationType } = req.body;
  logger.info(`Received /api/simulate request for SMILES: ${smiles}, Type: ${simulationType}`);

  if (!smiles) {
    return res.status(400).json({ message: 'SMILES string is required' });
  }
  if (!simulationType || !['binding', 'admet'].includes(simulationType)) {
    return res.status(400).json({ message: 'Invalid simulationType. Must be \'binding\' or \'admet\'.' });
  }

  try {
    // --- 1. Calculate Basic Properties --- 
    let properties = {};
    try {
      properties = await runRDKitScript(MOLECULAR_PROPERTIES_SCRIPT, [smiles]);
      logger.debug('Successfully calculated molecular properties');
    } catch (propError) {
      logger.warn(`Could not calculate properties for ${smiles}: ${propError.message}`);
      // Don't fail the whole simulation, just return empty properties
    }

    // --- 2. Perform Simulation based on Type --- 
    let simulationData = {};
    if (simulationType === 'admet') {
      try {
        simulationData = await runRDKitScript(ADMET_PREDICTION_SCRIPT, [smiles]);
        logger.debug('Successfully predicted ADMET properties');
      } catch (admetError) {
        logger.error(`ADMET prediction failed for ${smiles}: ${admetError.message}`);
        throw new Error(`ADMET prediction failed: ${admetError.message}`);
      }
    } else if (simulationType === 'binding') {
      // TODO: Implement actual binding affinity prediction/docking call
      // Requires target receptor info, etc.
      // Returning MOCK data for now
      logger.warn('Returning MOCK data for binding affinity simulation');
      simulationData = {
        'Dopamine Transporter': { score: Math.floor(Math.random() * 50) + 50, classification: 'Strong' }, // Random score 50-100
        'Norepinephrine Transporter': { score: Math.floor(Math.random() * 40) + 30, classification: 'Moderate' }, // Random score 30-70
        'Serotonin Transporter': { score: Math.floor(Math.random() * 30) + 1, classification: 'Weak' }, // Random score 1-30
        'D1 Receptor': { score: Math.floor(Math.random() * 25) + 1, classification: 'Weak' }, // Random score 1-25
        'D2 Receptor': { score: Math.floor(Math.random() * 40) + 20, classification: 'Moderate' } // Random score 20-60
      };
    }

    // --- 3. Format Response --- 
    // Attempt to get a molecule name (placeholder - needs implementation)
    // const moleculeName = getMoleculeName(smiles) || 'Simulated Molecule'; 
    const moleculeName = 'Simulated Molecule'; // Using placeholder for now

    const results = {
      type: simulationType,
      moleculeName: moleculeName, 
      smiles: smiles,
      properties: properties, // Include calculated properties
      ...(simulationType === 'binding' && { bindingAffinities: simulationData }),
      ...(simulationType === 'admet' && { admet: simulationData }),
    };

    // --- 4. Save and Respond --- 
    const resultId = uuidv4();
    const resultFile = path.join(resultsDir, `simulation_${resultId}.json`);
    fs.writeFileSync(resultFile, JSON.stringify({ id: resultId, timestamp: new Date().toISOString(), ...results }, null, 2));
    logger.info(`Simulation successful (ID: ${resultId}). Type: ${simulationType}.`);

    return res.json(results);

  } catch (error) {
    logger.error(`Error during /api/simulate (${simulationType}) for ${smiles}: ${error.message}`, { stack: error.stack });
    // Ensure we send a JSON response with a message field
    return res.status(500).json({
      message: `Error during ${simulationType} simulation: ${error.message}`,
    });
  }
});
// --- END NEW SIMULATE ENDPOINT ---

module.exports = router;
 