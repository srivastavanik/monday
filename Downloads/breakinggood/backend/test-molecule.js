// Script to test molecule generation and RDKit
const axios = require('axios');
const fs = require('fs');
const path = require('path');

// Fallback SMILES strings
const FALLBACK_SMILES = [
  'CN1C2=C(C=C(C=C2)Cl)N(C(=O)CC1=O)CC3=CC=C(C=C3)F',  // Flurazepam
  'C1=CC=C(C=C1)C2=COC3=CC(=CC(=C3)OC4=CC=CC=C4)C2=O',  // Flavanone
  'CC(C)(C)NCC(COC1=CC=CC2=CC=CC=C21)O',  // Propranolol
  'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  // Caffeine
  'CC(CS)C(=O)N1CCCC1C(=O)O'  // Captopril
];

// Test RDKit property calculation
async function testMolecules() {
  console.log('Testing molecule generation and RDKit integration');
  
  try {
    // Test each fallback molecule
    for (const smiles of FALLBACK_SMILES) {
      console.log(`Testing SMILES: ${smiles}`);
      
      try {
        // Call properties endpoint
        const propertiesResponse = await axios.post('http://localhost:5000/api/simulation/properties', {
          smiles: smiles
        });
        
        console.log('Properties response:', propertiesResponse.data);
        
        // Call ADMET endpoint
        const admetResponse = await axios.post('http://localhost:5000/api/simulation/admet', {
          smiles: smiles
        });
        
        console.log('ADMET response:', admetResponse.data);
      } catch (error) {
        console.error(`Error testing ${smiles}:`, error.message);
      }
    }
    
    // Now test the molecule generation process
    const generationResponse = await axios.post('http://localhost:5000/api/ai/generate-molecule', {
      requirements: 'Design a non-stimulant ADHD medication that targets dopamine transporters with minimal side effects.'
    });
    
    console.log('Generation response:', generationResponse.data);
    
    // Save the result for inspection
    fs.writeFileSync(path.join(__dirname, 'test-result.json'), 
      JSON.stringify(generationResponse.data, null, 2));
      
    console.log('Test complete! Check test-result.json for the output.');
  } catch (error) {
    console.error('Test failed:', error.message);
    if (error.response) {
      console.error('Response data:', error.response.data);
    }
  }
}

// Run the test
testMolecules();
