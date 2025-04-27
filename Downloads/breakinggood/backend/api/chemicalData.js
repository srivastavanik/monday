const express = require('express');
const axios = require('axios');
const router = express.Router();
const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');
const logger = require('../utils/logger');

// API Base URLs
const PUBCHEM_API_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
const CHEMBL_API_BASE = 'https://www.ebi.ac.uk/chembl/api/data';

// Python script to use chembl_webresource_client
const PYTHON_CHEMBL_SCRIPT = path.join(__dirname, '../scripts/chembl_search.py');

// Create scripts directory if it doesn't exist
const scriptsDir = path.join(__dirname, '../scripts');
if (!fs.existsSync(scriptsDir)) {
  fs.mkdirSync(scriptsDir);
}

// Create the Python script for ChEMBL integration
fs.writeFileSync(
  PYTHON_CHEMBL_SCRIPT,
  `
import sys
import json
import traceback
from chembl_webresource_client.new_client import new_client

def search_chembl(query_type, query_value, limit=10):
    try:
        results = []
        
        if query_type == 'molecule':
            molecule = new_client.molecule
            if query_value.startswith('CHEMBL'):
                res = molecule.filter(chembl_id=query_value).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
            else:
                res = molecule.filter(pref_name__icontains=query_value).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
            results = list(res)[:limit]
            
        elif query_type == 'target':
            target = new_client.target
            if query_value.startswith('CHEMBL'):
                res = target.filter(target_chembl_id=query_value)
            else:
                res = target.filter(pref_name__icontains=query_value)
            results = list(res)[:limit]
            
        elif query_type == 'activity':
            activity = new_client.activity
            if query_value.startswith('CHEMBL'):
                res = activity.filter(molecule_chembl_id=query_value).only(['molecule_chembl_id', 'target_chembl_id', 'standard_type', 'standard_value', 'standard_units'])
            else:
                res = activity.filter(pref_name__icontains=query_value).only(['molecule_chembl_id', 'target_chembl_id', 'standard_type', 'standard_value', 'standard_units'])
            results = list(res)[:limit]
            
        elif query_type == 'similarity':
            similarity = new_client.similarity
            res = similarity.filter(smiles=query_value, similarity=85).only(['molecule_chembl_id', 'pref_name', 'molecule_structures', 'similarity'])
            results = list(res)[:limit]
            
        elif query_type == 'substructure':
            substructure = new_client.substructure
            res = substructure.filter(smiles=query_value).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
            results = list(res)[:limit]
            
        return json.dumps(results)
    except Exception as e:
        return json.dumps({"error": str(e), "traceback": traceback.format_exc()})

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(json.dumps({"error": "Not enough arguments provided"}))
        sys.exit(1)
        
    query_type = sys.argv[1]
    query_value = sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    
    print(search_chembl(query_type, query_value, limit))
`
);

// Cache results in memory
const resultsCache = new Map();
const CACHE_TTL = 60 * 60 * 1000; // 1 hour in milliseconds

// Function to call the Python script with ChEMBL webresource client
const callChEMBLScript = (queryType, queryValue, limit = 10) => {
  return new Promise((resolve, reject) => {
    const cacheKey = `${queryType}:${queryValue}:${limit}`;
    
    // Check cache first
    if (resultsCache.has(cacheKey)) {
      const cachedData = resultsCache.get(cacheKey);
      if (Date.now() - cachedData.timestamp < CACHE_TTL) {
        return resolve(cachedData.data);
      }
      resultsCache.delete(cacheKey); // Expired cache entry
    }
    
    const pythonProcess = spawn('python', [
      PYTHON_CHEMBL_SCRIPT,
      queryType,
      queryValue,
      limit.toString()
    ]);
    
    let resultData = '';
    let errorData = '';
    
    pythonProcess.stdout.on('data', (data) => {
      resultData += data.toString();
    });
    
    pythonProcess.stderr.on('data', (data) => {
      errorData += data.toString();
    });
    
    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        logger.error(`ChEMBL script process exited with code ${code}: ${errorData}`);
        return reject(new Error(`ChEMBL query failed: ${errorData}`));
      }
      
      try {
        const result = JSON.parse(resultData);
        // Cache the result
        resultsCache.set(cacheKey, {
          data: result,
          timestamp: Date.now()
        });
        resolve(result);
      } catch (error) {
        logger.error(`Failed to parse ChEMBL result: ${error.message}`);
        reject(error);
      }
    });
  });
};

// Search compounds by name or similar structure
router.get('/search', async (req, res) => {
  try {
    const { query, source = 'pubchem', type = 'name', limit = 10 } = req.query;
    
    if (!query) {
      return res.status(400).json({ error: 'Search query is required' });
    }
    
    let results = [];
    
    if (source === 'pubchem') {
      // PubChem search
      if (type === 'name') {
        // Search by name
        const response = await axios.get(`${PUBCHEM_API_BASE}/compound/name/${encodeURIComponent(query)}/JSON`);
        
        if (response.data && response.data.PC_Compounds) {
          results = await processPublChemResults(response.data.PC_Compounds, limit);
        }
      } else if (type === 'smiles') {
        // Search by SMILES
        const response = await axios.get(`${PUBCHEM_API_BASE}/compound/smiles/${encodeURIComponent(query)}/JSON`);
        
        if (response.data && response.data.PC_Compounds) {
          results = await processPublChemResults(response.data.PC_Compounds, limit);
        }
      } else if (type === 'similarity') {
        // Search by similarity
        const response = await axios.get(
          `${PUBCHEM_API_BASE}/compound/similarity/${encodeURIComponent(query)}/JSON?Threshold=80&MaxRecords=${limit}`
        );
        
        if (response.data && response.data.PC_Compounds) {
          results = await processPublChemResults(response.data.PC_Compounds, limit);
        }
      }
    } else if (source === 'chembl') {
      // ChEMBL search
      if (type === 'name') {
        // Search by name
        const response = await axios.get(`${CHEMBL_API_BASE}/molecule?pref_name__icontains=${encodeURIComponent(query)}&limit=${limit}`);
        
        if (response.data && response.data.molecules) {
          results = processChemblResults(response.data.molecules, limit);
        }
      } else if (type === 'smiles') {
        // Search by structure
        const response = await axios.get(`${CHEMBL_API_BASE}/molecule?molecule_structures__canonical_smiles__flexmatch=${encodeURIComponent(query)}&limit=${limit}`);
        
        if (response.data && response.data.molecules) {
          results = processChemblResults(response.data.molecules, limit);
        }
      } else if (type === 'similarity') {
        // Search by similarity
        const response = await axios.get(
          `${CHEMBL_API_BASE}/similarity/${encodeURIComponent(query)}/80?limit=${limit}`
        );
        
        if (response.data && response.data.molecules) {
          results = processChemblResults(response.data.molecules, limit);
        }
      }
    }
    
    return res.json({ results });
    
  } catch (error) {
    console.error('Chemical search error:', error);
    return res.status(500).json({ 
      error: 'Error searching chemical database',
      details: error.message
    });
  }
});

// Get compound details by ID
router.get('/:source/:id', async (req, res) => {
  try {
    const { source, id } = req.params;
    
    if (!id) {
      return res.status(400).json({ error: 'Compound ID is required' });
    }
    
    let compoundDetails = null;
    
    if (source === 'pubchem') {
      const response = await axios.get(`${PUBCHEM_API_BASE}/compound/cid/${id}/JSON`);
      
      if (response.data && response.data.PC_Compounds && response.data.PC_Compounds.length > 0) {
        const compound = response.data.PC_Compounds[0];
        
        // Get properties in a separate call
        const propsResponse = await axios.get(
          `${PUBCHEM_API_BASE}/compound/cid/${id}/property/MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,ExactMass,MonoisotopicMass,Complexity,Charge/JSON`
        );
        
        let properties = {};
        if (propsResponse.data && propsResponse.data.PropertyTable && 
            propsResponse.data.PropertyTable.Properties && 
            propsResponse.data.PropertyTable.Properties.length > 0) {
          properties = propsResponse.data.PropertyTable.Properties[0];
        }
        
        // Get 2D structure
        const smiles = await getPublChemSmiles(id);
        
        compoundDetails = {
          id,
          source: 'pubchem',
          name: getPublChemValue(compound, 'IUPAC Name'),
          smiles,
          formula: properties.MolecularFormula || '',
          molecularWeight: properties.MolecularWeight || 0,
          logP: properties.XLogP || 0,
          tpsa: properties.TPSA || 0,
          hDonors: properties.HBondDonorCount || 0,
          hAcceptors: properties.HBondAcceptorCount || 0,
          rotatableBonds: properties.RotatableBondCount || 0,
          synonyms: getPublChemSynonyms(compound),
          imageUrl: `https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${id}`
        };
      }
    } else if (source === 'chembl') {
      const response = await axios.get(`${CHEMBL_API_BASE}/molecule/${id}`);
      
      if (response.data) {
        const molecule = response.data;
        const structures = molecule.molecule_structures || {};
        
        compoundDetails = {
          id,
          source: 'chembl',
          name: molecule.pref_name || 'Unknown',
          smiles: structures.canonical_smiles || '',
          formula: molecule.molecule_properties?.full_molformula || '',
          molecularWeight: molecule.molecule_properties?.full_mwt || 0,
          logP: molecule.molecule_properties?.alogp || 0,
          tpsa: molecule.molecule_properties?.psa || 0,
          hDonors: molecule.molecule_properties?.hbd || 0,
          hAcceptors: molecule.molecule_properties?.hba || 0,
          rotatableBonds: molecule.molecule_properties?.rtb || 0,
          synonyms: molecule.molecule_synonyms?.map(syn => syn.synonym) || [],
          imageUrl: `https://www.ebi.ac.uk/chembl/api/data/image/${id}.svg`
        };
      }
    }
    
    if (!compoundDetails) {
      return res.status(404).json({ error: 'Compound not found' });
    }
    
    return res.json(compoundDetails);
    
  } catch (error) {
    console.error('Error fetching compound details:', error);
    return res.status(500).json({ 
      error: 'Error fetching compound details',
      details: error.message
    });
  }
});

// Get bioactivity data for a compound
router.get('/:source/:id/bioactivity', async (req, res) => {
  try {
    const { source, id } = req.params;
    const { limit = 10 } = req.query;
    
    if (!id) {
      return res.status(400).json({ error: 'Compound ID is required' });
    }
    
    let bioactivityData = [];
    
    if (source === 'pubchem') {
      // For PubChem, we need to find ChEMBL ID first or use the PubChem assay data
      try {
        const response = await axios.get(`${PUBCHEM_API_BASE}/compound/cid/${id}/assaysummary/JSON`);
        
        if (response.data && response.data.AssaySummaries && response.data.AssaySummaries.length > 0) {
          bioactivityData = response.data.AssaySummaries
            .filter(assay => assay.ActiveOutcomeCount > 0)
            .slice(0, limit)
            .map(assay => ({
              assayId: assay.AID,
              assayName: assay.Name || 'Unknown',
              activeOutcomeCount: assay.ActiveOutcomeCount,
              totalOutcomeCount: assay.TotalOutcomeCount,
              source: 'PubChem BioAssay',
              url: `https://pubchem.ncbi.nlm.nih.gov/assay/${assay.AID}`
            }));
        }
      } catch (error) {
        console.error('Error fetching PubChem assay data:', error);
      }
    } else if (source === 'chembl') {
      try {
        const response = await axios.get(`${CHEMBL_API_BASE}/molecule/${id}/bioactivities`);
        
        if (response.data && response.data.bioactivities) {
          bioactivityData = response.data.bioactivities
            .slice(0, limit)
            .map(activity => ({
              assayId: activity.assay_chemblid,
              assayDescription: activity.assay_description || 'Unknown',
              targetName: activity.target_pref_name || 'Unknown',
              organism: activity.target_organism || 'Unknown',
              type: activity.standard_type || 'Unknown',
              relation: activity.standard_relation || '=',
              value: activity.standard_value || 0,
              units: activity.standard_units || '',
              source: 'ChEMBL',
              url: `https://www.ebi.ac.uk/chembl/assay_report_card/${activity.assay_chemblid}/`
            }));
        }
      } catch (error) {
        console.error('Error fetching ChEMBL bioactivity data:', error);
      }
    }
    
    return res.json({ bioactivityData });
    
  } catch (error) {
    console.error('Error fetching bioactivity data:', error);
    return res.status(500).json({ 
      error: 'Error fetching bioactivity data',
      details: error.message
    });
  }
});

// Helper functions
async function processPublChemResults(compounds, limit) {
  const results = [];
  
  for (let i = 0; i < Math.min(compounds.length, limit); i++) {
    const compound = compounds[i];
    const cid = compound.id.id.cid;
    
    // Get SMILES in a separate call
    const smiles = await getPublChemSmiles(cid);
    
    results.push({
      id: cid,
      source: 'pubchem',
      name: getPublChemValue(compound, 'IUPAC Name'),
      smiles,
      synonyms: getPublChemSynonyms(compound).slice(0, 3),
      imageUrl: `https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid=${cid}`
    });
  }
  
  return results;
}

function processChemblResults(molecules, limit) {
  return molecules
    .slice(0, limit)
    .map(molecule => {
      const structures = molecule.molecule_structures || {};
      
      return {
        id: molecule.molecule_chembl_id,
        source: 'chembl',
        name: molecule.pref_name || 'Unknown',
        smiles: structures.canonical_smiles || '',
        synonyms: (molecule.molecule_synonyms || [])
          .slice(0, 3)
          .map(syn => syn.synonym),
        imageUrl: `https://www.ebi.ac.uk/chembl/api/data/image/${molecule.molecule_chembl_id}.svg`
      };
    });
}

async function getPublChemSmiles(cid) {
  try {
    const response = await axios.get(`${PUBCHEM_API_BASE}/compound/cid/${cid}/property/CanonicalSMILES/JSON`);
    
    if (response.data && 
        response.data.PropertyTable && 
        response.data.PropertyTable.Properties && 
        response.data.PropertyTable.Properties.length > 0) {
      return response.data.PropertyTable.Properties[0].CanonicalSMILES || '';
    }
    
    return '';
  } catch (error) {
    console.error('Error fetching SMILES:', error);
    return '';
  }
}

function getPublChemValue(compound, key) {
  if (compound.props) {
    for (const prop of compound.props) {
      if (prop.urn.label === key) {
        if (prop.value.sval) {
          return prop.value.sval;
        } else if (prop.value.fval) {
          return prop.value.fval;
        } else if (prop.value.ival) {
          return prop.value.ival;
        }
      }
    }
  }
  
  return 'Unknown';
}

function getPublChemSynonyms(compound) {
  if (compound.props) {
    for (const prop of compound.props) {
      if (prop.urn.label === 'Synonym') {
        if (prop.value.stringValueList) {
          return prop.value.stringValueList.string || [];
        }
      }
    }
  }
  
  return [];
}

// Search ChEMBL database
router.get('/chembl', async (req, res) => {
  try {
    const { query, type = 'molecule', limit = 10 } = req.query;
    
    if (!query) {
      return res.status(400).json({ error: 'Query parameter is required' });
    }
    
    const results = await callChEMBLScript(type, query, limit);
    
    if (results.error) {
      logger.error(`ChEMBL API error: ${results.error}`);
      return res.status(500).json({ error: results.error });
    }
    
    res.json(results);
  } catch (error) {
    logger.error(`Error in ChEMBL search: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Get compound details from ChEMBL
router.get('/chembl/:chemblId', async (req, res) => {
  try {
    const { chemblId } = req.params;
    const results = await callChEMBLScript('molecule', chemblId);
    
    if (results.error) {
      return res.status(500).json({ error: results.error });
    }
    
    if (!results || results.length === 0) {
      return res.status(404).json({ error: 'Compound not found' });
    }
    
    res.json(results[0]);
  } catch (error) {
    logger.error(`Error getting ChEMBL compound: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Search by similarity
router.post('/similarity', async (req, res) => {
  try {
    const { smiles, threshold = 85 } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES structure is required' });
    }
    
    const results = await callChEMBLScript('similarity', smiles);
    res.json(results);
  } catch (error) {
    logger.error(`Error in similarity search: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Search by substructure
router.post('/substructure', async (req, res) => {
  try {
    const { smiles } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES structure is required' });
    }
    
    const results = await callChEMBLScript('substructure', smiles);
    res.json(results);
  } catch (error) {
    logger.error(`Error in substructure search: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Get bioactivity data
router.get('/bioactivity/:chemblId', async (req, res) => {
  try {
    const { chemblId } = req.params;
    const results = await callChEMBLScript('activity', chemblId);
    
    if (results.error) {
      return res.status(500).json({ error: results.error });
    }
    
    res.json(results);
  } catch (error) {
    logger.error(`Error getting bioactivity data: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

// Get drug-target interactions
router.get('/drug-target/:chemblId', async (req, res) => {
  try {
    const { chemblId } = req.params;
    const activityResults = await callChEMBLScript('activity', chemblId);
    
    if (activityResults.error) {
      return res.status(500).json({ error: activityResults.error });
    }
    
    // Extract target IDs
    const targetIds = [...new Set(activityResults.map(item => item.target_chembl_id))];
    
    // Get target details
    const targetDetails = await Promise.all(
      targetIds.map(targetId => callChEMBLScript('target', targetId))
    );
    
    const interactions = targetIds.map((targetId, index) => ({
      targetId,
      targetDetails: targetDetails[index],
      activities: activityResults.filter(item => item.target_chembl_id === targetId)
    }));
    
    res.json(interactions);
  } catch (error) {
    logger.error(`Error getting drug-target interactions: ${error.message}`);
    res.status(500).json({ error: error.message });
  }
});

module.exports = router; 