const express = require('express');
const router = express.Router();
const axios = require('axios'); // Required for calling other internal APIs
const logger = require('../utils/logger');

// Generate full regulatory report using other services and Claude
router.post('/report', async (req, res) => {
  try {
    const { 
      smiles,
      drugClass = 'CNS stimulant',
      targetIndication = 'ADHD',
      primaryMechanism = 'dopamine/norepinephrine reuptake inhibition',
      novelMechanism = false,
      orphanDrug = false,
      fastTrack = false
    } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES string is required' });
    }
    
    logger.info(`Generating regulatory report for SMILES: ${smiles}`);

    // 1. Get Molecular Properties from Simulation API
    let properties = null;
    try {
      const propResponse = await axios.post('http://localhost:5000/api/simulation/properties', { smiles });
      properties = propResponse.data;
      if (properties.error) throw new Error(properties.error); // Handle errors from the properties call
      logger.info('Successfully fetched molecular properties.');
    } catch (error) {
      logger.error(`Error fetching properties for ${smiles}: ${error.message}`);
      return res.status(500).json({ error: 'Failed to calculate molecule properties', details: error.message });
    }

    // 2. Get ADMET Prediction from Simulation API
    let admet = null;
    try {
      const admetResponse = await axios.post('http://localhost:5000/api/simulation/admet', { smiles });
      admet = admetResponse.data;
      if (admet.error) throw new Error(admet.error);
      logger.info('Successfully fetched ADMET predictions.');
    } catch (error) {
      logger.error(`Error fetching ADMET for ${smiles}: ${error.message}`);
      // Continue without ADMET if it fails, Claude can still provide some analysis
      admet = { error: `Failed to predict ADMET: ${error.message}` };
    }

    // 3. Get Similar Approved Drugs (Example using ChEMBL similarity - could be more sophisticated)
    let similarDrugs = [];
    try {
      // Use the dedicated similarity search endpoint
      const similarityResponse = await axios.post('http://localhost:5000/api/similarity/search', { 
          query: smiles, 
          targets: [], // Ideally, search against an approved drug database target list
          threshold: 85, 
          limit: 5 
      }); 
      // This endpoint needs adjustment or a different target source to find *approved* drugs.
      // For now, we log the limitation and return an empty array or results from general ChEMBL search.
      logger.warn('Similar drug search currently uses general ChEMBL similarity, not specific approved drug database.');
      // If the similarity endpoint returns results:
      // similarDrugs = similarityResponse.data?.filter(d => d.pref_name).slice(0, 5) || []; 
      logger.info(`Found ${similarDrugs.length} potentially similar ChEMBL compounds.`);
    } catch (error) {
        logger.error(`Error fetching similar drugs: ${error.message}`);
        // Continue without similar drugs if search fails
    }

    // 4. Ask Claude to Generate the Regulatory Report based on gathered data
    const systemPrompt = `You are an expert in pharmaceutical regulatory affairs, medicinal chemistry, and drug development. Generate a comprehensive regulatory analysis report for the candidate molecule described below. Structure the report clearly.`;
    
    const userPrompt = `Please generate a regulatory analysis report for the following drug candidate:

**Molecule SMILES:** ${smiles}

**Calculated Properties:**
${JSON.stringify(properties, null, 2)}

**Predicted ADMET Profile:**
${JSON.stringify(admet, null, 2)}

**Target Indication:** ${targetIndication}
**Drug Class:** ${drugClass}
**Primary Mechanism:** ${primaryMechanism}
**Is Mechanism Novel?** ${novelMechanism}
**Potential Orphan Drug Status?** ${orphanDrug}
**Potential Fast Track Status?** ${fastTrack}

**Potentially Similar Approved/Known Drugs (for context):**
${similarDrugs.length > 0 ? similarDrugs.map(d => `- ${d.pref_name} (ChEMBL ID: ${d.molecule_chembl_id}, Similarity: ${d.similarity?.toFixed(2)})`).join('\n') : '- None found in basic similarity search.'}

Please include the following sections in your report:
1.  **Executive Summary:** Brief overview of the candidate, potential pathway, and key considerations.
2.  **Regulatory Pathway Estimation:** Estimated timelines (Preclinical, Phase I, II, III, FDA Review, Total), potential for special designations (Orphan, Fast Track, etc.), key milestones.
3.  **Key Risks and Challenges:** Based on properties, ADMET, mechanism, and drug class (e.g., safety concerns, CMC issues, scheduling).
4.  **Required Studies Overview:** Summary of typical studies needed for IND and NDA submission for this type of drug and indication.
5.  **Competitive Landscape Snapshot:** Briefly mention key existing treatments for ${targetIndication} and potential advantages/disadvantages of this candidate based on the provided data.
6.  **Production & CMC Considerations:** Brief comment on potential synthetic complexity challenges based on the structure (e.g., complex ring systems, stereocenters - inferred).

Provide realistic estimations and justifications based *only* on the provided data. Acknowledge limitations where data is missing.`;

    let analysisReportText = '';
    try {
      // Call Claude via the AI service endpoint
      const claudeResponse = await axios.post('http://localhost:5000/api/ai/ask', {
        question: userPrompt, // Send the full structured prompt
        context: `Regulatory analysis for SMILES: ${smiles}` // Provide context
      });
      analysisReportText = claudeResponse.data.response;
      logger.info('Successfully generated regulatory analysis report via Claude.');
    } catch (error) {
        logger.error(`Error calling Claude for regulatory analysis: ${error.message}`);
        analysisReportText = `Error: Could not generate AI analysis. ${error.message}`;
    }

    // Return the structured results including the AI report
    return res.json({
      smiles,
      parameters: { drugClass, targetIndication, primaryMechanism, novelMechanism, orphanDrug, fastTrack },
      properties,
      admet,
      similarDrugs, // Include the found similar drugs
      regulatoryReport: analysisReportText // The AI-generated text report
    });
    
  } catch (error) {
    logger.error(`Regulatory report generation error: ${error.message}`, { stack: error.stack });
    return res.status(500).json({ 
      error: 'Error generating regulatory report',
      details: error.message
    });
  }
});

// Other endpoints remain as placeholders or can be removed
router.post('/pathway', async (req, res) => {
   logger.warn('Endpoint /api/regulatory/pathway called - DEPRECATED. Use /api/regulatory/report instead.');
   res.status(404).json({ error: 'Endpoint deprecated. Use /api/regulatory/report instead.' });
});

router.post('/production-feasibility', async (req, res) => {
   logger.warn('Endpoint /api/regulatory/production-feasibility called - DEPRECATED. Use /api/regulatory/report instead.');
   res.status(404).json({ error: 'Endpoint deprecated. Use /api/regulatory/report instead.' });
});

router.post('/market-analysis', async (req, res) => {
   logger.warn('Endpoint /api/regulatory/market-analysis called - DEPRECATED. Use /api/regulatory/report instead.');
   res.status(404).json({ error: 'Endpoint deprecated. Use /api/regulatory/report instead.' });
});

module.exports = router; 