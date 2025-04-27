// Handler functions for MoleculeDesigner component

// Handler for selecting a molecule from ThinkingProcess component
const handleSelectFromThinking = async (smiles, simulationAPI, aiGeneratedMolecules, setAiGeneratedMolecules, setSelectedMolecule, aiRequestId, showSnackbar) => {
  console.log('Selected molecule from thinking process:', smiles);
  
  // If smiles is valid, create a molecule object and select it
  if (smiles && typeof smiles === 'string') {
    try {
      // Get molecule properties using RDKit
      const propertiesResponse = await simulationAPI.runRDKit({
        smiles: smiles,
        operation: 'descriptors'
      });
      
      const newMolecule = {
        id: Date.now(),
        name: `AI-Compound-${aiGeneratedMolecules.length + 1}`,
        smiles: smiles,
        properties: {
          molecularWeight: `${propertiesResponse.data.molecular_weight.toFixed(2)} g/mol`,
          logP: propertiesResponse.data.logp.toFixed(2),
          hBondDonors: propertiesResponse.data.num_h_donors,
          hBondAcceptors: propertiesResponse.data.num_h_acceptors,
          rotableBonds: propertiesResponse.data.num_rotatable_bonds,
          psa: `${propertiesResponse.data.tpsa.toFixed(1)} Å²`,
          formula: propertiesResponse.data.formula,
        },
        timestamp: new Date().toISOString(),
        aiGenerated: true,
        requestId: aiRequestId,
      };
      
      setAiGeneratedMolecules(prev => [...prev, newMolecule]);
      setSelectedMolecule(newMolecule);
      showSnackbar('Molecule selected from AI thinking process');
    } catch (err) {
      console.error('Error creating molecule from SMILES:', err);
      
      // Create basic molecule without properties
      const basicMolecule = {
        id: Date.now(),
        name: `AI-Compound-${aiGeneratedMolecules.length + 1}`,
        smiles: smiles,
        properties: { formula: 'Unknown' },
        timestamp: new Date().toISOString(),
        aiGenerated: true,
        requestId: aiRequestId,
      };
      
      setAiGeneratedMolecules(prev => [...prev, basicMolecule]);
      setSelectedMolecule(basicMolecule);
      showSnackbar('Molecule selected (with limited properties)');
    }
  } else {
    console.error('Invalid SMILES from ThinkingProcess:', smiles);
    showSnackbar('Failed to select molecule: Invalid structure', 'error');
  }
};

// Handler for saving a molecule from ThinkingProcess component
const handleSaveFromThinking = (molecule, savedMolecules, setSavedMolecules, showSnackbar) => {
  console.log('Saving molecule from thinking process:', molecule);
  
  if (!molecule || !molecule.smiles) {
    showSnackbar('Cannot save molecule: Invalid data', 'error');
    return;
  }
  
  // Check if molecule already exists in savedMolecules
  const exists = savedMolecules.some(m => 
    m.smiles === molecule.smiles || 
    (m.id && molecule.id && m.id === molecule.id)
  );
  
  if (!exists) {
    // Create a proper molecule object if just a SMILES string was provided
    const moleculeToSave = typeof molecule === 'string' ? 
      {
        id: Date.now(),
        name: `Saved-${savedMolecules.length + 1}`,
        smiles: molecule,
        properties: { formula: 'Unknown' },
        timestamp: new Date().toISOString(),
      } : molecule;
    
    const updatedSavedMolecules = [...savedMolecules, moleculeToSave];
    setSavedMolecules(updatedSavedMolecules);
    
    // Also save to localStorage
    try {
      localStorage.setItem('savedMolecules', JSON.stringify(updatedSavedMolecules));
    } catch (err) {
      console.error('Error saving to localStorage:', err);
    }
    
    showSnackbar(`${moleculeToSave.name || 'Molecule'} saved successfully`);
  } else {
    showSnackbar('This molecule is already saved', 'warning');
  }
};

export { handleSelectFromThinking, handleSaveFromThinking };
