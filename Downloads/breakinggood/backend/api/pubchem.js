const express = require('express');
const axios = require('axios');
const fs = require('fs');
const path = require('path');
const router = express.Router();
const logger = require('../utils/logger');

// PubChem API base URL
const PUBCHEM_BASE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';

// Create storage directory for chemical data
const chemicalDataDir = path.join(__dirname, '../data/chemical_data');
if (!fs.existsSync(chemicalDataDir)) {
  fs.mkdirSync(chemicalDataDir, { recursive: true });
}

// Get compound data by CID, name, or SMILES
router.get('/compound', async (req, res) => {
  try {
    const { identifier, identifierType = 'name', properties = 'MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey' } = req.query;
    
    if (!identifier) {
      return res.status(400).json({ error: 'Compound identifier is required' });
    }
    
    // Validate identifier type
    const validIdentifierTypes = ['cid', 'name', 'smiles', 'inchi', 'inchikey'];
    if (!validIdentifierTypes.includes(identifierType)) {
      return res.status(400).json({ error: `Identifier type must be one of: ${validIdentifierTypes.join(', ')}` });
    }
    
    // Build the PubChem API URL
    let requestUrl = `${PUBCHEM_BASE_URL}/compound/${identifierType}/${encodeURIComponent(identifier)}/property/${properties}/JSON`;
    
    // Special handling for SMILES
    if (identifierType === 'smiles') {
      requestUrl = `${PUBCHEM_BASE_URL}/compound/smiles/${encodeURIComponent(identifier)}/property/${properties}/JSON`;
    }
    
    // Query PubChem API
    const response = await axios.get(requestUrl);
    
    // Save the retrieved data
    const sanitizedIdentifier = identifier.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 50);
    const fileName = `pubchem_${identifierType}_${sanitizedIdentifier}.json`;
    const filePath = path.join(chemicalDataDir, fileName);
    
    fs.writeFileSync(filePath, JSON.stringify(response.data, null, 2));
    
    return res.json({
      source: 'PubChem',
      identifier,
      identifierType,
      data: response.data,
      savedTo: fileName
    });
    
  } catch (error) {
    logger.error(`Error fetching PubChem data: ${error.message}`);
    return res.status(500).json({
      error: 'Error fetching compound data from PubChem',
      details: error.message
    });
  }
});

// Get compound 2D image
router.get('/compound/image', async (req, res) => {
  try {
    const { identifier, identifierType = 'name', size = 300 } = req.query;
    
    if (!identifier) {
      return res.status(400).json({ error: 'Compound identifier is required' });
    }
    
    // Validate identifier type
    const validIdentifierTypes = ['cid', 'name', 'smiles', 'inchi', 'inchikey'];
    if (!validIdentifierTypes.includes(identifierType)) {
      return res.status(400).json({ error: `Identifier type must be one of: ${validIdentifierTypes.join(', ')}` });
    }
    
    // Build the PubChem API URL for image
    let requestUrl = `${PUBCHEM_BASE_URL}/compound/${identifierType}/${encodeURIComponent(identifier)}/PNG?image_size=${size}x${size}`;
    
    // Special handling for SMILES
    if (identifierType === 'smiles') {
      requestUrl = `${PUBCHEM_BASE_URL}/compound/smiles/${encodeURIComponent(identifier)}/PNG?image_size=${size}x${size}`;
    }
    
    // Query PubChem API and get image
    const response = await axios.get(requestUrl, { responseType: 'arraybuffer' });
    
    // Set appropriate headers and return the image
    res.set('Content-Type', 'image/png');
    return res.send(Buffer.from(response.data, 'binary'));
    
  } catch (error) {
    logger.error(`Error fetching compound image: ${error.message}`);
    return res.status(500).json({
      error: 'Error fetching compound image from PubChem',
      details: error.message
    });
  }
});

// Search compounds by query
router.get('/search', async (req, res) => {
  try {
    const { query, maxResults = 10 } = req.query;
    
    if (!query) {
      return res.status(400).json({ error: 'Search query is required' });
    }
    
    // Build the PubChem search URL
    const searchUrl = `${PUBCHEM_BASE_URL}/compound/name/${encodeURIComponent(query)}/cids/JSON`;
    
    // Query PubChem search API
    const searchResponse = await axios.get(searchUrl);
    
    if (!searchResponse.data.IdentifierList || !searchResponse.data.IdentifierList.CID || searchResponse.data.IdentifierList.CID.length === 0) {
      return res.json({
        query,
        results: []
      });
    }
    
    // Limit the number of results
    const cids = searchResponse.data.IdentifierList.CID.slice(0, maxResults);
    
    // Get detailed info for each CID
    const properties = 'MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey';
    const detailsUrl = `${PUBCHEM_BASE_URL}/compound/cid/${cids.join(',')}/property/${properties}/JSON`;
    
    const detailsResponse = await axios.get(detailsUrl);
    
    return res.json({
      query,
      totalResults: searchResponse.data.IdentifierList.CID.length,
      results: detailsResponse.data.PropertyTable.Properties
    });
    
  } catch (error) {
    logger.error(`Error searching PubChem: ${error.message}`);
    return res.status(500).json({
      error: 'Error searching compounds on PubChem',
      details: error.message
    });
  }
});

// Get compound assay data (bioactivity)
router.get('/bioactivity', async (req, res) => {
  try {
    const { cid, target } = req.query;
    
    if (!cid) {
      return res.status(400).json({ error: 'Compound CID is required' });
    }
    
    // Build the PubChem assay URL
    let assayUrl;
    if (target) {
      // Search for target-specific bioactivity
      assayUrl = `${PUBCHEM_BASE_URL}/compound/cid/${cid}/assaysummary/JSON?target=${encodeURIComponent(target)}`;
    } else {
      // Get all bioactivity data
      assayUrl = `${PUBCHEM_BASE_URL}/compound/cid/${cid}/assaysummary/JSON`;
    }
    
    // Query PubChem assay API
    const assayResponse = await axios.get(assayUrl);
    
    // Save the retrieved data
    const fileName = `pubchem_bioactivity_${cid}${target ? '_' + target.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 30) : ''}.json`;
    const filePath = path.join(chemicalDataDir, fileName);
    
    fs.writeFileSync(filePath, JSON.stringify(assayResponse.data, null, 2));
    
    return res.json({
      source: 'PubChem',
      cid,
      target: target || null,
      data: assayResponse.data,
      savedTo: fileName
    });
    
  } catch (error) {
    logger.error(`Error fetching bioactivity data: ${error.message}`);
    return res.status(500).json({
      error: 'Error fetching bioactivity data from PubChem',
      details: error.message
    });
  }
});

module.exports = router;
