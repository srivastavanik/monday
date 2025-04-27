// Fix for Structure Editor page

// Debounce function to limit how often a function is called
function debounce(func, wait) {
  let timeout;
  return function() {
    const context = this;
    const args = arguments;
    clearTimeout(timeout);
    timeout = setTimeout(() => {
      func.apply(context, args);
    }, wait);
  };
}

document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying structure editor fixes...');
  
  // Fix structure editor
  fixStructureEditor();
  setInterval(fixStructureEditor, 2000);
});

// Main structure editor fix function
function fixStructureEditor() {
  // Find the structure editor tab
  const structureEditorTab = document.querySelector('[class*="structureEditor"], #structure-editor');
  if (!structureEditorTab) return;
  
  // Find SMILES input field
  const smilesInput = structureEditorTab.querySelector('input, textarea');
  if (!smilesInput) return;
  
  // Find molecule preview container
  let previewContainer = structureEditorTab.querySelector('.molecule-preview-container');
  if (!previewContainer) {
    const existingContainer = structureEditorTab.querySelector('[class*="preview"]');
    
    if (existingContainer && existingContainer.textContent.includes('preview')) {
      // Create a proper container
      previewContainer = document.createElement('div');
      previewContainer.className = 'molecule-preview-container';
      previewContainer.style.cssText = 'border:1px solid #e0e0e0;border-radius:4px;height:300px;margin-top:16px;';
      
      // Replace existing container with safety checks
      if (existingContainer && existingContainer.parentNode) {
        try {
          existingContainer.parentNode.replaceChild(previewContainer, existingContainer);
        } catch (error) {
          console.warn('Error replacing container:', error);
          // Alternative approach - insert after instead of replacing
          if (existingContainer.parentNode) {
            existingContainer.parentNode.insertBefore(previewContainer, existingContainer.nextSibling);
            // Hide the original container instead of removing it
            existingContainer.style.display = 'none';
          }
        }
      } else {
        // If we can't find a proper place, just append to the tab
        structureEditorTab.appendChild(previewContainer);
      }
    }
  }
  
  // Setup visualization of molecule based on SMILES input
  if (smilesInput && previewContainer) {
    // Only add listener if not already added
    if (!smilesInput.dataset.hasListener) {
      smilesInput.dataset.hasListener = 'true';
      
      // Listen for input changes
      smilesInput.addEventListener('input', debounce(function(e) {
        const smiles = e.target.value.trim();
        if (smiles.length > 5) {
          // Visualize the molecule
          previewContainer.innerHTML = '';
          createMoleculeVisualization(previewContainer, smiles);
        }
      }, 500));
      
      // Find update/import buttons
      const updateButton = structureEditorTab.querySelector('button:not(.MuiButton-contained)');
      if (updateButton) {
        updateButton.addEventListener('click', function() {
          const smiles = smilesInput.value.trim();
          if (smiles.length > 5) {
            previewContainer.innerHTML = '';
            createMoleculeVisualization(previewContainer, smiles);
          }
        });
      }
      
      // Handle sample molecule
      const sampleSmiles = 'CC1=C(C(=O)OC1=O)C2=CC=CC=C2';
      if (!smilesInput.value) {
        smilesInput.value = sampleSmiles;
        previewContainer.innerHTML = '';
        createMoleculeVisualization(previewContainer, sampleSmiles);
      }
    }
  }
  
  // Fix the molecule name and save button
  const nameInput = structureEditorTab.querySelector('input[type="text"]:not([class*="smiles"])');
  const saveButton = structureEditorTab.querySelector('button[class*="save"], button:contains("SAVE")');
  
  if (nameInput && saveButton && !saveButton.dataset.hasListener) {
    saveButton.dataset.hasListener = 'true';
    
    saveButton.addEventListener('click', function() {
      const smiles = smilesInput.value.trim();
      const name = nameInput.value.trim() || `Molecule-${Date.now()}`;
      
      if (smiles.length > 5) {
        saveMolecule(name, smiles);
      }
    });
  }
}

// Create a molecule visualization (this should match the function in molecule-visualization-fix.js)
function createMoleculeVisualization(container, smilesText) {
  try {
    // Extract SMILES if formatted with label
    let smiles = smilesText;
    if (smilesText.includes('SMILES:')) {
      const match = smilesText.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
      if (match) smiles = match[1].trim();
    }
    
    // Create a unique ID
    const viewerId = 'mol-viewer-' + Math.random().toString(36).substring(2, 9);
    
    // Create a container for the viewer
    container.innerHTML = `
      <div id="${viewerId}" style="width:100%;height:100%;position:relative;">
        <div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);">
          <div class="spinner" style="border:4px solid rgba(0,0,0,0.1);border-radius:50%;border-top:4px solid #3498db;width:30px;height:30px;animation:spin 1s linear infinite;"></div>
        </div>
      </div>
      <style>
        @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
      </style>
    `;
    
    // Get the viewer element
    const viewerElement = document.getElementById(viewerId);
    
    // Make API call to get 3D structure
    fetch('/api/simulation/3d-structure', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: smiles })
    })
    .then(response => response.json())
    .then(data => {
      if (!data.molblock) throw new Error('No molblock received');
      
      // Create viewer
      const viewer = $3Dmol.createViewer(viewerElement, {
        backgroundColor: 'white'
      });
      
      // Add the model
      viewer.addModel(data.molblock, 'mol');
      
      // Style the molecule
      viewer.setStyle({}, { stick: {} });
      viewer.addSurface($3Dmol.SurfaceType.VDW, {
        opacity: 0.6,
        color: 'lightblue'
      });
      
      // Center and render
      viewer.zoomTo();
      viewer.render();
      
      // Store for cleanup
      viewerElement.viewer = viewer;
    })
    .catch(error => {
      console.error('Error visualizing molecule:', error);
      // Add a null check before setting innerHTML
      if (viewerElement) {
        viewerElement.innerHTML = `<div style="color:red;padding:10px;text-align:center;">Error: ${error.message}</div>`;
      } else {
        console.warn('Viewer element not found when handling error');
        // Try to find the container and update it instead
        if (container) {
          container.innerHTML = `<div style="color:red;padding:10px;text-align:center;">Error: ${error.message}</div>`;
        }
      }
    });
  } catch (error) {
    console.error('Error creating visualization:', error);
    container.innerHTML = `<div style="color:red;padding:10px;text-align:center;">Error creating molecule visualization</div>`;
  }
}

// Save a molecule to localStorage
function saveMolecule(name, smiles) {
  try {
    // Create molecule object
    const molecule = {
      id: Date.now(),
      name: name,
      smiles: smiles,
      timestamp: new Date().toISOString()
    };
    
    // Get existing saved molecules
    let savedMolecules = [];
    try {
      const saved = localStorage.getItem('savedMolecules');
      if (saved) {
        savedMolecules = JSON.parse(saved);
      }
    } catch (err) {
      console.error('Error loading saved molecules:', err);
    }
    
    // Add new molecule
    savedMolecules.push(molecule);
    
    // Save back to localStorage
    localStorage.setItem('savedMolecules', JSON.stringify(savedMolecules));
    
    // Show success message
    showToast(`Saved molecule: ${name}`);
    
    // Update saved molecules tab if it's visible
    updateSavedMoleculesTab();
  } catch (error) {
    console.error('Error saving molecule:', error);
    showToast('Error saving molecule', 'error');
  }
}

// Update the saved molecules tab
function updateSavedMoleculesTab() {
  // Find the saved molecules tab
  const savedTab = document.querySelector('[class*="savedMolecules"], #saved-molecules');
  if (!savedTab) return;
  
  try {
    // Get saved molecules
    const saved = localStorage.getItem('savedMolecules');
    if (!saved) return;
    
    const molecules = JSON.parse(saved);
    if (!molecules || molecules.length === 0) return;
    
    // Create a container for the molecules
    const container = document.createElement('div');
    container.style.cssText = 'display:flex;flex-wrap:wrap;gap:16px;margin-top:16px;';
    
    // Create a card for each molecule
    molecules.forEach(molecule => {
      const card = document.createElement('div');
      card.className = 'saved-molecule-card';
      card.style.cssText = 'width:300px;border:1px solid #e0e0e0;border-radius:4px;overflow:hidden;';
      
      card.innerHTML = `
        <div style="padding:12px;background-color:#f5f5f5;">
          <h3 style="margin:0;">${molecule.name}</h3>
          <div style="font-size:0.8em;color:#666;">${new Date(molecule.timestamp).toLocaleString()}</div>
        </div>
        <div class="molecule-preview" style="height:200px;"></div>
        <div style="padding:12px;">
          <pre style="overflow:auto;background:#f9f9f9;padding:8px;font-size:12px;">${molecule.smiles}</pre>
          <div style="display:flex;gap:8px;margin-top:12px;">
            <button class="select-btn" style="flex:1;padding:8px;background:#2196f3;color:white;border:none;border-radius:4px;cursor:pointer;">SELECT</button>
            <button class="delete-btn" style="flex:1;padding:8px;background:#f44336;color:white;border:none;border-radius:4px;cursor:pointer;">DELETE</button>
          </div>
        </div>
      `;
      
      container.appendChild(card);
    });
    
    // Replace existing content
    savedTab.innerHTML = '<h2>Saved Molecules</h2>';
    savedTab.appendChild(container);
    
    // Add visualization to each molecule
    const previewContainers = container.querySelectorAll('.molecule-preview');
    previewContainers.forEach((previewContainer, index) => {
      if (index < molecules.length) {
        createMoleculeVisualization(previewContainer, molecules[index].smiles);
      }
    });
    
    // Add event listeners to buttons
    const selectButtons = container.querySelectorAll('.select-btn');
    selectButtons.forEach((button, index) => {
      button.addEventListener('click', function() {
        if (index < molecules.length) {
          // Find the structure editor tab and update it
          const structureEditorTab = document.querySelector('[class*="structureEditor"], #structure-editor');
          if (structureEditorTab) {
            const smilesInput = structureEditorTab.querySelector('input, textarea');
            const nameInput = structureEditorTab.querySelector('input[type="text"]:not([class*="smiles"])');
            
            if (smilesInput) smilesInput.value = molecules[index].smiles;
            if (nameInput) nameInput.value = molecules[index].name;
            
            // Switch to structure editor tab
            const tabButtons = document.querySelectorAll('[role="tab"]');
            tabButtons.forEach(tabButton => {
              if (tabButton.textContent.includes('STRUCTURE') || tabButton.textContent.includes('Editor')) {
                tabButton.click();
              }
            });
          }
        }
      });
    });
    
    const deleteButtons = container.querySelectorAll('.delete-btn');
    deleteButtons.forEach((button, index) => {
      button.addEventListener('click', function() {
        if (index < molecules.length) {
          // Remove from saved molecules
          molecules.splice(index, 1);
          localStorage.setItem('savedMolecules', JSON.stringify(molecules));
          
          // Update the view
          updateSavedMoleculesTab();
          showToast('Molecule deleted');
        }
      });
    });
  } catch (error) {
    console.error('Error updating saved molecules tab:', error);
  }
}

// Show a toast message
function showToast(message, type = 'success') {
  // Create toast container if it doesn't exist
  let toastContainer = document.querySelector('.toast-container');
  if (!toastContainer) {
    toastContainer = document.createElement('div');
    toastContainer.className = 'toast-container';
    toastContainer.style.cssText = 'position:fixed;bottom:20px;right:20px;z-index:9999;';
    document.body.appendChild(toastContainer);
  }
  
  // Create toast
  const toast = document.createElement('div');
  toast.className = `toast ${type}`;
  toast.style.cssText = `
    margin-top:10px;
    padding:12px 16px;
    border-radius:4px;
    color:white;
    font-size:14px;
    box-shadow:0 2px 5px rgba(0,0,0,0.2);
    animation:fadeIn 0.3s, fadeOut 0.3s 2.7s;
    background-color:${type === 'success' ? '#4caf50' : '#f44336'};
  `;
  
  toast.innerHTML = message;
  
  // Add to container
  toastContainer.appendChild(toast);
  
  // Remove after 3 seconds
  setTimeout(() => {
    toast.remove();
  }, 3000);
}

// Debounce function to prevent too many calls
function debounce(func, wait) {
  let timeout;
  return function(...args) {
    const context = this;
    clearTimeout(timeout);
    timeout = setTimeout(() => func.apply(context, args), wait);
  };
}
