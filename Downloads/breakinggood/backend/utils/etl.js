const fs = require('fs');
const path = require('path');
const axios = require('axios');
const { exec } = require('child_process');
const logger = require('./logger');

// Data directories
const dataDir = path.join(__dirname, '../data');
const articlesDir = path.join(dataDir, 'articles');
const chemicalDataDir = path.join(dataDir, 'chemical_data');
const etlResultsDir = path.join(dataDir, 'etl_results');

// Ensure directories exist
for (const dir of [dataDir, articlesDir, chemicalDataDir, etlResultsDir]) {
  if (!fs.existsSync(dir)) {
    fs.mkdirSync(dir, { recursive: true });
  }
}

/**
 * Extract literature data from PubMed and PMC
 * @param {Object} options - Options for extraction
 * @param {string} options.query - Search query
 * @param {number} options.maxResults - Maximum number of results to extract
 * @param {boolean} options.includeBioC - Whether to include BioC format
 * @param {boolean} options.includeFullText - Whether to include full text
 * @returns {Promise<Object>} Extraction results
 */
const extractLiterature = async (options) => {
  try {
    const { query, maxResults = 100, includeBioC = true, includeFullText = true } = options;
    
    // Search PubMed for articles matching the query
    const searchUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi';
    const searchResponse = await axios.get(searchUrl, {
      params: {
        db: 'pubmed',
        term: query,
        retmax: maxResults,
        retmode: 'json',
        sort: 'relevance'
      }
    });
    
    const idList = searchResponse.data.esearchresult.idlist || [];
    
    if (idList.length === 0) {
      return { success: false, error: 'No articles found matching the query' };
    }
    
    // Fetch article details
    const fetchUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
    const fetchResponse = await axios.get(fetchUrl, {
      params: {
        db: 'pubmed',
        id: idList.join(','),
        retmode: 'xml'
      }
    });
    
    // Write raw XML data
    const xmlFilePath = path.join(etlResultsDir, `pubmed_${new Date().toISOString().replace(/[:.]/g, '-')}.xml`);
    fs.writeFileSync(xmlFilePath, fetchResponse.data);
    
    // Extract PMC IDs for full text retrieval
    const pmcIdRegex = /PMC\d+/g;
    const pmcIds = fetchResponse.data.match(pmcIdRegex) || [];
    
    // Convert PMC IDs to the format expected by the API (remove 'PMC' prefix)
    const formattedPmcIds = pmcIds
      .map(id => id.replace('PMC', ''))
      .filter((value, index, self) => self.indexOf(value) === index); // Remove duplicates
    
    // Results to return
    const results = {
      query,
      totalFound: idList.length,
      articlesProcessed: idList.length,
      xmlSavedTo: xmlFilePath,
      pmcIds: formattedPmcIds,
      fullTextArticles: [],
      biocArticles: []
    };
    
    // Fetch full text if requested
    if (includeFullText && formattedPmcIds.length > 0) {
      for (const pmcid of formattedPmcIds) {
        try {
          // Fetch full text in XML format
          const fullTextUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
          const fullTextResponse = await axios.get(fullTextUrl, {
            params: {
              db: 'pmc',
              id: pmcid,
              retmode: 'xml'
            }
          });
          
          const fullTextPath = path.join(articlesDir, `PMC${pmcid}_fulltext.xml`);
          fs.writeFileSync(fullTextPath, fullTextResponse.data);
          
          results.fullTextArticles.push({
            pmcid,
            savedTo: fullTextPath
          });
        } catch (ftError) {
          logger.error(`Error fetching full text for PMC${pmcid}: ${ftError.message}`);
        }
      }
    }
    
    // Fetch BioC format if requested
    if (includeBioC && formattedPmcIds.length > 0) {
      for (const pmcid of formattedPmcIds) {
        try {
          // Fetch article in BioC format
          const biocUrl = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi';
          const biocResponse = await axios.get(biocUrl, {
            params: {
              datatype: 'pmc',
              format: 'json',
              id: pmcid
            }
          });
          
          const biocPath = path.join(articlesDir, `PMC${pmcid}_bioc.json`);
          fs.writeFileSync(biocPath, JSON.stringify(biocResponse.data, null, 2));
          
          results.biocArticles.push({
            pmcid,
            savedTo: biocPath
          });
          
          // Extract and save annotations
          extractAnnotations(biocResponse.data, pmcid);
        } catch (biocError) {
          logger.error(`Error fetching BioC format for PMC${pmcid}: ${biocError.message}`);
        }
      }
    }
    
    // Save extraction results
    const resultsSummaryPath = path.join(etlResultsDir, `literature_extraction_${new Date().toISOString().replace(/[:.]/g, '-')}.json`);
    fs.writeFileSync(resultsSummaryPath, JSON.stringify(results, null, 2));
    
    return { success: true, results };
  } catch (error) {
    logger.error(`Literature extraction error: ${error.message}`);
    return { success: false, error: error.message };
  }
};

/**
 * Extract annotations from BioC format data
 * @param {Object} biocData - BioC format data
 * @param {string} pmcid - PMC ID
 * @returns {Object} Extraction results
 */
const extractAnnotations = (biocData, pmcid) => {
  try {
    const annotations = {
      receptors: [],
      chemicals: [],
      genes: [],
      diseases: [],
      methods: [],
      outcomes: []
    };
    
    if (biocData.documents) {
      for (const document of biocData.documents) {
        if (document.passages) {
          for (const passage of document.passages) {
            // Extract section type (if available)
            const sectionType = passage.infons && passage.infons.section_type ? passage.infons.section_type : 'unknown';
            
            // Process annotations if available
            if (passage.annotations) {
              for (const annotation of passage.annotations) {
                const type = annotation.infons && annotation.infons.type ? annotation.infons.type.toLowerCase() : 'unknown';
                const text = annotation.text;
                
                // Categorize annotations
                if (type === 'chemical' || type === 'compound') {
                  annotations.chemicals.push({ type, text, section: sectionType });
                } else if (type === 'gene' || type === 'protein') {
                  annotations.genes.push({ type, text, section: sectionType });
                } else if (type === 'disease' || type === 'disorder') {
                  annotations.diseases.push({ type, text, section: sectionType });
                } else if (type === 'method' || type === 'technique') {
                  annotations.methods.push({ type, text, section: sectionType });
                }
                
                // Look for receptor keywords
                if (text.toLowerCase().includes('receptor') || 
                    text.toLowerCase().includes('transporters')) {
                  annotations.receptors.push({ type, text, section: sectionType });
                }
                
                // Look for outcome measure keywords
                if (sectionType === 'results' || sectionType === 'conclusions') {
                  if (text.toLowerCase().includes('outcome') || 
                      text.toLowerCase().includes('efficacy') || 
                      text.toLowerCase().includes('improvement')) {
                    annotations.outcomes.push({ type, text, section: sectionType });
                  }
                }
              }
            }
            
            // Text mining for receptor targets, methods, and outcomes
            // even if not formally annotated
            const passageText = passage.text.toLowerCase();
            
            // Check for receptor targets
            const receptorKeywords = ['receptor', 'transporter', 'dopamine', 'norepinephrine', 
                                    'serotonin', 'gaba', 'glutamate'];
            
            for (const keyword of receptorKeywords) {
              if (passageText.includes(keyword)) {
                // Find the sentence containing the keyword
                const sentences = passage.text.split(/\.s+/);
                for (const sentence of sentences) {
                  if (sentence.toLowerCase().includes(keyword)) {
                    annotations.receptors.push({
                      type: 'text_mining',
                      text: sentence.trim(),
                      section: sectionType,
                      keyword
                    });
                    break; // Just get the first sentence with this keyword
                  }
                }
              }
            }
            
            // Check for experimental methods
            const methodKeywords = ['assay', 'analysis', 'experiment', 'measurement', 
                                  'technique', 'procedure', 'test', 'model'];
            
            if (sectionType === 'methods' || sectionType === 'materials') {
              for (const keyword of methodKeywords) {
                if (passageText.includes(keyword)) {
                  // Find the sentence containing the keyword
                  const sentences = passage.text.split(/\.s+/);
                  for (const sentence of sentences) {
                    if (sentence.toLowerCase().includes(keyword)) {
                      annotations.methods.push({
                        type: 'text_mining',
                        text: sentence.trim(),
                        section: sectionType,
                        keyword
                      });
                      break; // Just get the first sentence with this keyword
                    }
                  }
                }
              }
            }
            
            // Check for outcome measures
            const outcomeKeywords = ['outcome', 'efficacy', 'improvement', 'reduction', 
                                    'increase', 'decrease', 'significant'];
            
            if (sectionType === 'results' || sectionType === 'conclusions') {
              for (const keyword of outcomeKeywords) {
                if (passageText.includes(keyword)) {
                  // Find the sentence containing the keyword
                  const sentences = passage.text.split(/\.s+/);
                  for (const sentence of sentences) {
                    if (sentence.toLowerCase().includes(keyword)) {
                      annotations.outcomes.push({
                        type: 'text_mining',
                        text: sentence.trim(),
                        section: sectionType,
                        keyword
                      });
                      break; // Just get the first sentence with this keyword
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    // Remove duplicates
    for (const category in annotations) {
      annotations[category] = annotations[category].filter((item, index, self) => 
        index === self.findIndex(t => t.text === item.text)
      );
    }
    
    // Save annotations
    const annotationsPath = path.join(articlesDir, `PMC${pmcid}_annotations.json`);
    fs.writeFileSync(annotationsPath, JSON.stringify(annotations, null, 2));
    
    return annotations;
  } catch (error) {
    logger.error(`Annotation extraction error for PMC${pmcid}: ${error.message}`);
    return {};
  }
};

/**
 * Extract chemical data from PubChem and ChEMBL
 * @param {Object} options - Options for extraction
 * @param {string[]} options.smiles - SMILES strings to extract data for
 * @param {string[]} options.names - Chemical names to extract data for
 * @param {boolean} options.includePubChem - Whether to include PubChem data
 * @param {boolean} options.includeChEMBL - Whether to include ChEMBL data
 * @returns {Promise<Object>} Extraction results
 */
const extractChemicalData = async (options) => {
  try {
    const { smiles = [], names = [], includePubChem = true, includeChEMBL = true } = options;
    
    if (smiles.length === 0 && names.length === 0) {
      return { success: false, error: 'No SMILES or names provided' };
    }
    
    const results = {
      totalProcessed: smiles.length + names.length,
      pubchemResults: [],
      chemblResults: []
    };
    
    // Process SMILES strings
    for (const smilesStr of smiles) {
      try {
        // Extract from PubChem
        if (includePubChem) {
          const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smilesStr)}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey/JSON`;
          const pubchemResponse = await axios.get(pubchemUrl);
          
          const fileName = `pubchem_smiles_${smilesStr.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 30)}.json`;
          const filePath = path.join(chemicalDataDir, fileName);
          
          fs.writeFileSync(filePath, JSON.stringify(pubchemResponse.data, null, 2));
          
          results.pubchemResults.push({
            smiles: smilesStr,
            savedTo: fileName,
            data: pubchemResponse.data.PropertyTable.Properties[0]
          });
        }
        
        // Extract from ChEMBL
        if (includeChEMBL) {
          const chemblUrl = `https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__canonical_smiles=${encodeURIComponent(smilesStr)}`;
          const chemblResponse = await axios.get(chemblUrl);
          
          if (chemblResponse.data.molecules && chemblResponse.data.molecules.length > 0) {
            const fileName = `chembl_smiles_${smilesStr.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 30)}.json`;
            const filePath = path.join(chemicalDataDir, fileName);
            
            fs.writeFileSync(filePath, JSON.stringify(chemblResponse.data, null, 2));
            
            // Get bioactivity data if available
            let bioactivityData = null;
            const firstMolecule = chemblResponse.data.molecules[0];
            
            if (firstMolecule && firstMolecule.molecule_chembl_id) {
              try {
                const bioactivityUrl = `https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=${firstMolecule.molecule_chembl_id}`;
                const bioactivityResponse = await axios.get(bioactivityUrl);
                
                const bioactivityFileName = `chembl_bioactivity_${firstMolecule.molecule_chembl_id}.json`;
                const bioactivityFilePath = path.join(chemicalDataDir, bioactivityFileName);
                
                fs.writeFileSync(bioactivityFilePath, JSON.stringify(bioactivityResponse.data, null, 2));
                
                bioactivityData = {
                  savedTo: bioactivityFileName,
                  count: bioactivityResponse.data.activities ? bioactivityResponse.data.activities.length : 0
                };
              } catch (bioErr) {
                logger.error(`Error fetching bioactivity for ${firstMolecule.molecule_chembl_id}: ${bioErr.message}`);
              }
            }
            
            results.chemblResults.push({
              smiles: smilesStr,
              savedTo: fileName,
              chemblId: firstMolecule.molecule_chembl_id,
              bioactivity: bioactivityData
            });
          }
        }
      } catch (smilesErr) {
        logger.error(`Error processing SMILES ${smilesStr}: ${smilesErr.message}`);
      }
    }
    
    // Process chemical names
    for (const name of names) {
      try {
        // Extract from PubChem
        if (includePubChem) {
          const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey/JSON`;
          const pubchemResponse = await axios.get(pubchemUrl);
          
          const fileName = `pubchem_name_${name.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 30)}.json`;
          const filePath = path.join(chemicalDataDir, fileName);
          
          fs.writeFileSync(filePath, JSON.stringify(pubchemResponse.data, null, 2));
          
          results.pubchemResults.push({
            name,
            savedTo: fileName,
            data: pubchemResponse.data.PropertyTable.Properties[0]
          });
        }
        
        // Extract from ChEMBL
        if (includeChEMBL) {
          const chemblUrl = `https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name=${encodeURIComponent(name)}`;
          const chemblResponse = await axios.get(chemblUrl);
          
          if (chemblResponse.data.molecules && chemblResponse.data.molecules.length > 0) {
            const fileName = `chembl_name_${name.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 30)}.json`;
            const filePath = path.join(chemicalDataDir, fileName);
            
            fs.writeFileSync(filePath, JSON.stringify(chemblResponse.data, null, 2));
            
            // Get bioactivity data if available
            let bioactivityData = null;
            const firstMolecule = chemblResponse.data.molecules[0];
            
            if (firstMolecule && firstMolecule.molecule_chembl_id) {
              try {
                const bioactivityUrl = `https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=${firstMolecule.molecule_chembl_id}`;
                const bioactivityResponse = await axios.get(bioactivityUrl);
                
                const bioactivityFileName = `chembl_bioactivity_${firstMolecule.molecule_chembl_id}.json`;
                const bioactivityFilePath = path.join(chemicalDataDir, bioactivityFileName);
                
                fs.writeFileSync(bioactivityFilePath, JSON.stringify(bioactivityResponse.data, null, 2));
                
                bioactivityData = {
                  savedTo: bioactivityFileName,
                  count: bioactivityResponse.data.activities ? bioactivityResponse.data.activities.length : 0
                };
              } catch (bioErr) {
                logger.error(`Error fetching bioactivity for ${firstMolecule.molecule_chembl_id}: ${bioErr.message}`);
              }
            }
            
            results.chemblResults.push({
              name,
              savedTo: fileName,
              chemblId: firstMolecule.molecule_chembl_id,
              bioactivity: bioactivityData
            });
          }
        }
      } catch (nameErr) {
        logger.error(`Error processing name ${name}: ${nameErr.message}`);
      }
    }
    
    // Save extraction results
    const resultsSummaryPath = path.join(etlResultsDir, `chemical_extraction_${new Date().toISOString().replace(/[:.]/g, '-')}.json`);
    fs.writeFileSync(resultsSummaryPath, JSON.stringify(results, null, 2));
    
    return { success: true, results };
  } catch (error) {
    logger.error(`Chemical data extraction error: ${error.message}`);
    return { success: false, error: error.message };
  }
};

/**
 * Run a complete ETL pipeline for literature and chemical data
 * @param {Object} options - ETL options
 * @returns {Promise<Object>} ETL results
 */
const runEtlPipeline = async (options) => {
  try {
    const {
      literatureQuery,
      maxLiteratureResults = 100,
      includeBioC = true,
      includeFullText = true,
      smiles = [],
      chemicalNames = [],
      includePubChem = true,
      includeChEMBL = true
    } = options;
    
    const results = {
      timestamp: new Date().toISOString(),
      literature: null,
      chemicalData: null
    };
    
    // Extract literature data if query provided
    if (literatureQuery) {
      logger.info(`Starting literature extraction for query: ${literatureQuery}`);
      
      const literatureResults = await extractLiterature({
        query: literatureQuery,
        maxResults: maxLiteratureResults,
        includeBioC,
        includeFullText
      });
      
      results.literature = literatureResults;
      
      // If literature extraction found chemicals, add them to chemical extraction
      if (literatureResults.success && literatureResults.results.biocArticles.length > 0) {
        // Extract chemical names from annotations
        const extractedChemicals = new Set();
        
        for (const article of literatureResults.results.biocArticles) {
          try {
            const annotationsPath = path.join(articlesDir, `PMC${article.pmcid}_annotations.json`);
            
            if (fs.existsSync(annotationsPath)) {
              const annotations = JSON.parse(fs.readFileSync(annotationsPath, 'utf8'));
              
              if (annotations.chemicals && annotations.chemicals.length > 0) {
                for (const chemical of annotations.chemicals) {
                  extractedChemicals.add(chemical.text);
                }
              }
            }
          } catch (annoErr) {
            logger.error(`Error extracting chemicals from annotations: ${annoErr.message}`);
          }
        }
        
        // Add extracted chemicals to the list
        if (extractedChemicals.size > 0) {
          chemicalNames.push(...Array.from(extractedChemicals));
        }
      }
    }
    
    // Extract chemical data if SMILES or names provided
    if (smiles.length > 0 || chemicalNames.length > 0) {
      logger.info(`Starting chemical data extraction for ${smiles.length} SMILES and ${chemicalNames.length} names`);
      
      const chemicalResults = await extractChemicalData({
        smiles,
        names: chemicalNames,
        includePubChem,
        includeChEMBL
      });
      
      results.chemicalData = chemicalResults;
    }
    
    // Save the complete ETL results
    const etlResultsPath = path.join(etlResultsDir, `etl_pipeline_${new Date().toISOString().replace(/[:.]/g, '-')}.json`);
    fs.writeFileSync(etlResultsPath, JSON.stringify(results, null, 2));
    
    return { success: true, results };
  } catch (error) {
    logger.error(`ETL pipeline error: ${error.message}`);
    return { success: false, error: error.message };
  }
};

module.exports = {
  extractLiterature,
  extractAnnotations,
  extractChemicalData,
  runEtlPipeline
};
