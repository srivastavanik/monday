const express = require('express');
const axios = require('axios');
const fs = require('fs');
const path = require('path');
const router = express.Router();
const logger = require('../utils/logger');

// PMC Open Access BioC API base URL
const PMC_BIOC_BASE_URL = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi';
const OA_SERVICE_BASE_URL = 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi';

// Create storage directory for downloaded articles
const articlesDir = path.join(__dirname, '../data/articles');
if (!fs.existsSync(articlesDir)) {
  fs.mkdirSync(articlesDir, { recursive: true });
}

// Get article in BioC format by PMC ID
router.get('/pmc/:pmcid', async (req, res) => {
  try {
    const { pmcid } = req.params;
    const { format = 'json' } = req.query;
    
    if (!pmcid) {
      return res.status(400).json({ error: 'PMC ID is required' });
    }
    
    // Format validation
    if (!['xml', 'json'].includes(format)) {
      return res.status(400).json({ error: 'Format must be xml or json' });
    }
    
    // Query PMC BioC API
    const response = await axios.get(PMC_BIOC_BASE_URL, {
      params: {
        datatype: 'pmc',
        format: format,
        id: pmcid
      },
      responseType: format === 'xml' ? 'text' : 'json'
    });
    
    // Store the article data
    const fileName = `PMC${pmcid}.${format}`;
    const filePath = path.join(articlesDir, fileName);
    
    if (format === 'xml') {
      fs.writeFileSync(filePath, response.data);
    } else {
      fs.writeFileSync(filePath, JSON.stringify(response.data, null, 2));
    }
    
    return res.json({
      pmcid,
      format,
      data: response.data,
      savedTo: fileName
    });
    
  } catch (error) {
    logger.error(`Error fetching BioC article: ${error.message}`);
    return res.status(500).json({
      error: 'Error fetching article in BioC format',
      details: error.message
    });
  }
});

// Query OA Web Service for updated articles
router.get('/oa/updated', async (req, res) => {
  try {
    const { from, format = 'json' } = req.query;
    
    if (!from) {
      return res.status(400).json({ error: 'From date is required (YYYY-MM-DD)' });
    }
    
    // Query OA Web Service API
    const response = await axios.get(OA_SERVICE_BASE_URL, {
      params: {
        format: 'json',
        from: from
      }
    });
    
    return res.json(response.data);
    
  } catch (error) {
    logger.error(`Error querying OA Web Service: ${error.message}`);
    return res.status(500).json({
      error: 'Error querying OA Web Service',
      details: error.message
    });
  }
});

// Download resources from OA Web Service
router.get('/oa/download', async (req, res) => {
  try {
    const { id, resource = 'pdf' } = req.query;
    
    if (!id) {
      return res.status(400).json({ error: 'Article ID is required' });
    }
    
    if (!['pdf', 'tgz'].includes(resource)) {
      return res.status(400).json({ error: 'Resource must be pdf or tgz' });
    }
    
    // Query OA Web Service for the resource URL
    const response = await axios.get(OA_SERVICE_BASE_URL, {
      params: {
        id: id,
        format: 'json'
      }
    });
    
    // Extract resource URL
    let resourceUrl = null;
    if (resource === 'pdf') {
      resourceUrl = response.data.records[0]?.pdfUrl;
    } else {
      resourceUrl = response.data.records[0]?.tgzUrl;
    }
    
    if (!resourceUrl) {
      return res.status(404).json({ error: `No ${resource} resource found for this article` });
    }
    
    // Download the resource
    const resourceResponse = await axios.get(resourceUrl, { responseType: 'stream' });
    
    // Save the resource
    const fileName = `${id}.${resource}`;
    const filePath = path.join(articlesDir, fileName);
    const writer = fs.createWriteStream(filePath);
    
    resourceResponse.data.pipe(writer);
    
    return new Promise((resolve, reject) => {
      writer.on('finish', () => {
        resolve(res.json({
          id,
          resource,
          savedTo: fileName
        }));
      });
      writer.on('error', err => {
        reject(res.status(500).json({
          error: 'Error saving resource',
          details: err.message
        }));
      });
    });
    
  } catch (error) {
    logger.error(`Error downloading resource: ${error.message}`);
    return res.status(500).json({
      error: 'Error downloading resource',
      details: error.message
    });
  }
});

// Extract and save annotations from BioC articles
router.post('/extract-annotations', async (req, res) => {
  try {
    const { pmcid, annotationTypes = ['chemical', 'disease', 'gene'] } = req.body;
    
    if (!pmcid) {
      return res.status(400).json({ error: 'PMC ID is required' });
    }
    
    // Check if article exists, if not, fetch it
    const fileName = `PMC${pmcid}.json`;
    const filePath = path.join(articlesDir, fileName);
    
    let articleData;
    
    if (!fs.existsSync(filePath)) {
      // Fetch the article first
      const response = await axios.get(PMC_BIOC_BASE_URL, {
        params: {
          datatype: 'pmc',
          format: 'json',
          id: pmcid
        }
      });
      
      articleData = response.data;
      fs.writeFileSync(filePath, JSON.stringify(articleData, null, 2));
    } else {
      articleData = JSON.parse(fs.readFileSync(filePath, 'utf8'));
    }
    
    // Extract annotations
    const annotations = [];
    
    if (articleData.documents) {
      for (const document of articleData.documents) {
        for (const passage of document.passages) {
          if (passage.annotations) {
            for (const annotation of passage.annotations) {
              if (annotationTypes.includes(annotation.infons.type)) {
                annotations.push({
                  type: annotation.infons.type,
                  text: annotation.text,
                  identifier: annotation.infons.identifier || null,
                  location: {
                    offset: annotation.locations[0].offset,
                    length: annotation.locations[0].length
                  },
                  passage: passage.text
                });
              }
            }
          }
        }
      }
    }
    
    // Save annotations
    const annotationsFileName = `PMC${pmcid}_annotations.json`;
    const annotationsFilePath = path.join(articlesDir, annotationsFileName);
    fs.writeFileSync(annotationsFilePath, JSON.stringify(annotations, null, 2));
    
    return res.json({
      pmcid,
      annotationTypes,
      totalAnnotations: annotations.length,
      annotations,
      savedTo: annotationsFileName
    });
    
  } catch (error) {
    logger.error(`Error extracting annotations: ${error.message}`);
    return res.status(500).json({
      error: 'Error extracting annotations',
      details: error.message
    });
  }
});

module.exports = router;
