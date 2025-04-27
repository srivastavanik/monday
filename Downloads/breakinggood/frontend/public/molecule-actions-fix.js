// Fix for molecule select/save functionality across tabs

document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying molecule actions fixes...');
  
  // Apply fixes to select/save buttons
  fixMoleculeActions();
  setInterval(fixMoleculeActions, 1500);
});

// Main function to fix molecule actions
function fixMoleculeActions() {
  // Find all molecule cards across the application
  findAndFixMoleculeCards();
  
  // Set up global event listeners for molecule selection sharing
  setupGlobalEventListeners();
  
  // Replace AI with Claude
  replaceAIWithClaude();
}

// Find and fix all molecule cards
function findAndFixMoleculeCards() {
  // Find all molecule cards across the application
  const moleculeCards = document.querySelectorAll('.molecule-card, [class*="moleculeCard"], [class*="MoleculeCard"]');
  const generatedCards = document.querySelectorAll('[class*="Generated"] .MuiPaper-root, [class*="generated"] .MuiPaper-root');
  const allCards = [...moleculeCards, ...generatedCards];
  
  allCards.forEach(card => {
    // Skip if already processed
    if (card.dataset.actionsFixed === 'true') return;
    
    // Find SMILES within the card
    const smiles = findSmilesInCard(card);
    if (!smiles) return;
    
    // Find existing or add select/save buttons
    fixSelectSaveButtons(card, smiles);
    
    // Mark as processed
    card.dataset.actionsFixed = 'true';
  });
}

// Find SMILES in a molecule card
function findSmilesInCard(card) {
  // Look for SMILES in code or pre elements
  let smilesElem = card.querySelector('code, pre, [class*="code"], [class*="smiles"]');
  if (smilesElem && smilesElem.textContent && isSmilesFormat(smilesElem.textContent)) {
    return extractSmiles(smilesElem.textContent);
  }
  
  // Look in attributes
  if (card.dataset.smiles && isSmilesFormat(card.dataset.smiles)) {
    return card.dataset.smiles;
  }
  
  // Check data attributes of child elements
  const children = card.querySelectorAll('*');
  for (const child of children) {
    if (child.dataset.smiles && isSmilesFormat(child.dataset.smiles)) {
      return child.dataset.smiles;
    }
  }
  
  // Search text content of descendants
  const potentialElements = card.querySelectorAll('p, div, span');
  for (const elem of potentialElements) {
    if (elem.textContent && isSmilesFormat(elem.textContent)) {
      return extractSmiles(elem.textContent);
    }
  }
  
  return null;
}

// Fix or add select/save buttons to a card
function fixSelectSaveButtons(card, smiles) {
  // Check if card already has buttons
  let selectButton = card.querySelector('[class*="select"], button:nth-of-type(1)');
  let saveButton = card.querySelector('[class*="save"], button:nth-of-type(2)');
  let buttonsContainer = card.querySelector('[class*="actions"], [class*="buttons"]');
  
  // If no buttons container found, create one
  if (!buttonsContainer) {
    buttonsContainer = document.createElement('div');
    buttonsContainer.className = 'molecule-actions';
    buttonsContainer.style.cssText = 'display:flex;gap:8px;margin-top:12px;padding:0 12px 12px;';
    
    // Find a good place to insert
    const insertTarget = card.querySelector('[class*="content"], [class*="CardContent"]') || card;
    insertTarget.appendChild(buttonsContainer);
  }
  
  // Create select button if not found
  if (!selectButton) {
    selectButton = document.createElement('button');
    selectButton.className = 'select-molecule-btn';
    selectButton.textContent = 'SELECT';
    selectButton.style.cssText = 'flex:1;padding:8px;background:#2196f3;color:white;border:none;border-radius:4px;cursor:pointer;font-weight:500;';
    buttonsContainer.appendChild(selectButton);
  }
  
  // Create save button if not found
  if (!saveButton) {
    saveButton = document.createElement('button');
    saveButton.className = 'save-molecule-btn';
    saveButton.textContent = 'SAVE';
    saveButton.style.cssText = 'flex:1;padding:8px;background:#4caf50;color:white;border:none;border-radius:4px;cursor:pointer;font-weight:500;';
    buttonsContainer.appendChild(saveButton);
  }
  
  // Ensure buttons have proper event listeners
  if (!selectButton.dataset.hasListener) {
    selectButton.dataset.hasListener = 'true';
    selectButton.addEventListener('click', function(event) {
      event.preventDefault();
      event.stopPropagation();
      selectMolecule(smiles, card);
    });
  }
  
  if (!saveButton.dataset.hasListener) {
    saveButton.dataset.hasListener = 'true';
    saveButton.addEventListener('click', function(event) {
      event.preventDefault();
      event.stopPropagation();
      saveMolecule(smiles, card);
    });
  }
}

// Select a molecule across the application
function selectMolecule(smiles, sourceCard) {
  console.log('Selecting molecule:', smiles.substring(0, 20) + '...');
  
  try {
    // Highlight the selected card
    highlightSelectedCard(sourceCard);
    
    // Store in global state for other parts of the application
    window.selectedMolecule = {
      smiles: smiles,
      timestamp: new Date().toISOString()
    };
    
    // Dispatch custom event for other components
    const selectEvent = new CustomEvent('molecule:selected', {
      detail: { smiles: smiles, sourceCard: sourceCard }
    });
    document.dispatchEvent(selectEvent);
    
    // Update relevant UI elements
    updateMoleculeDesigner(smiles);
    updateStructureEditor(smiles);
    updateMoleculeViewer(smiles);
    
    // Show feedback
    showToast('Molecule selected');
  } catch (error) {
    console.error('Error selecting molecule:', error);
    showToast('Error selecting molecule', 'error');
  }
}

// Highlight the selected card
function highlightSelectedCard(card) {
  // Remove highlight from other cards
  const allCards = document.querySelectorAll('.molecule-card, [class*="moleculeCard"], [class*="MoleculeCard"], [class*="Generated"] .MuiPaper-root');
  allCards.forEach(c => {
    c.style.borderColor = '';
    c.style.boxShadow = '';
  });
  
  // Highlight this card
  card.style.borderColor = '#2196f3';
  card.style.boxShadow = '0 0 8px rgba(33, 150, 243, 0.5)';
}

// Save a molecule
function saveMolecule(smiles, sourceCard) {
  console.log('Saving molecule:', smiles.substring(0, 20) + '...');
  
  try {
    // Generate name from card or create a default
    let name = 'Molecule-' + Math.floor(Math.random() * 10000);
    
    // Try to find a better name from the card
    const nameElem = sourceCard.querySelector('h3, h4, [class*="title"]');
    if (nameElem && nameElem.textContent) {
      name = nameElem.textContent.trim();
    }
    
    // Create molecule object
    const molecule = {
      id: Date.now(),
      name: name,
      smiles: smiles,
      timestamp: new Date().toISOString()
    };
    
    // Store in localStorage
    const savedMolecules = JSON.parse(localStorage.getItem('savedMolecules') || '[]');
    savedMolecules.push(molecule);
    localStorage.setItem('savedMolecules', JSON.stringify(savedMolecules));
    
    // Dispatch custom event
    const saveEvent = new CustomEvent('molecule:saved', {
      detail: { molecule: molecule, sourceCard: sourceCard }
    });
    document.dispatchEvent(saveEvent);
    
    // Update UI
    updateSavedMoleculesTab();
    
    // Show feedback
    showToast(`Molecule saved as "${name}"`);
  } catch (error) {
    console.error('Error saving molecule:', error);
    showToast('Error saving molecule', 'error');
  }
}

// Set up global event listeners
function setupGlobalEventListeners() {
  // Only set up once
  if (window.moleculeActionsListenersInitialized) return;
  window.moleculeActionsListenersInitialized = true;
  
  // Listen for tab changes to update relevant UIs
  const tabs = document.querySelectorAll('[role="tab"]');
  tabs.forEach(tab => {
    if (!tab.dataset.hasListener) {
      tab.dataset.hasListener = 'true';
      tab.addEventListener('click', function() {
        setTimeout(() => {
          if (window.selectedMolecule) {
            updateMoleculeViewer(window.selectedMolecule.smiles);
            updateStructureEditor(window.selectedMolecule.smiles);
          }
          updateSavedMoleculesTab();
        }, 100);
      });
    }
  });
}

// Update molecule designer with selected molecule
function updateMoleculeDesigner(smiles) {
  // Find the molecule designer component
  const designer = document.querySelector('[class*="MoleculeDesigner"], [class*="designer"]');
  if (!designer) return;
  
  // Find SMILES input or relevant field
  const inputs = designer.querySelectorAll('input, textarea');
  let smilesInput = null;
  
  for (const input of inputs) {
    // Look for clues this is a SMILES input
    if (input.id && input.id.toLowerCase().includes('smiles')) {
      smilesInput = input;
      break;
    }
    if (input.placeholder && input.placeholder.toLowerCase().includes('smiles')) {
      smilesInput = input;
      break;
    }
    if (input.name && input.name.toLowerCase().includes('smiles')) {
      smilesInput = input;
      break;
    }
    // Last resort
    if (input.type === 'text' || input.tagName === 'TEXTAREA') {
      smilesInput = input;
    }
  }
  
  if (smilesInput) {
    // Update input
    smilesInput.value = smiles;
    
    // Manually dispatch input event to trigger any listeners
    const event = new Event('input', { bubbles: true });
    smilesInput.dispatchEvent(event);
  }
}

// Update structure editor with selected molecule
function updateStructureEditor(smiles) {
  // Find the structure editor tab
  const editorTab = document.querySelector('[class*="structureEditor"], #structure-editor');
  if (!editorTab) return;
  
  // Find SMILES input
  const smilesInput = editorTab.querySelector('input, textarea');
  if (smilesInput) {
    // Update input
    smilesInput.value = smiles;
    
    // Manually dispatch input event to trigger visualization
    const event = new Event('input', { bubbles: true });
    smilesInput.dispatchEvent(event);
    
    // Also try to update preview
    const previewContainer = editorTab.querySelector('.molecule-preview-container, [class*="preview"]');
    if (previewContainer) {
      setTimeout(() => {
        // This function should be defined in structure-editor-fix.js
        if (typeof createMoleculeVisualization === 'function') {
          createMoleculeVisualization(previewContainer, smiles);
        }
      }, 100);
    }
  }
}

// Update molecule viewer with selected molecule
function updateMoleculeViewer(smiles) {
  // Find main molecule viewer
  const viewers = document.querySelectorAll('.MoleculeViewer3D, [class*="moleculeViewer"], [class*="Viewer3D"]');
  
  viewers.forEach(viewer => {
    // Skip if inside a card that might be showing a different molecule
    const isInCard = viewer.closest('.molecule-card, [class*="moleculeCard"], [class*="MoleculeCard"]');
    if (isInCard) return;
    
    // Create visualization
    if (typeof createMoleculeVisualization === 'function') {
      createMoleculeVisualization(viewer, smiles);
    }
  });
}

// Update saved molecules tab
function updateSavedMoleculesTab() {
  try {
    // Find the saved molecules tab
    const savedTab = document.querySelector('[class*="savedMolecules"], #saved-molecules');
    if (!savedTab) return;
    
    // Get saved molecules
    const savedData = localStorage.getItem('savedMolecules');
    if (!savedData) return;
    
    const molecules = JSON.parse(savedData);
    if (!molecules || molecules.length === 0) return;
    
    // Create a container for the saved molecules
    let container = savedTab.querySelector('.saved-molecules-container');
    if (!container) {
      container = document.createElement('div');
      container.className = 'saved-molecules-container';
      container.style.cssText = 'display:flex;flex-wrap:wrap;gap:16px;margin-top:16px;';
      
      // Add title if not already present
      let title = savedTab.querySelector('h2, h3');
      if (!title) {
        title = document.createElement('h2');
        title.textContent = 'Saved Molecules';
        savedTab.appendChild(title);
      }
      
      savedTab.appendChild(container);
    }
    
    // Clear container and add molecule cards
    container.innerHTML = '';
    
    molecules.forEach(molecule => {
      const card = document.createElement('div');
      card.className = 'molecule-card';
      card.style.cssText = 'width:300px;border:1px solid #e0e0e0;border-radius:4px;overflow:hidden;';
      
      card.innerHTML = `
        <div style="padding:12px;background-color:#f5f5f5;">
          <h3 style="margin:0;">${molecule.name}</h3>
          <div style="font-size:0.8em;color:#666;">${new Date(molecule.timestamp).toLocaleString()}</div>
        </div>
        <div class="molecule-preview" style="height:200px;"></div>
        <div style="padding:12px;">
          <pre style="overflow:auto;background:#f9f9f9;padding:8px;font-size:12px;">${molecule.smiles}</pre>
        </div>
      `;
      
      // Add buttons container
      const buttonsContainer = document.createElement('div');
      buttonsContainer.style.cssText = 'display:flex;gap:8px;margin-top:12px;padding:0 12px 12px;';
      
      // Add select button
      const selectButton = document.createElement('button');
      selectButton.textContent = 'SELECT';
      selectButton.style.cssText = 'flex:1;padding:8px;background:#2196f3;color:white;border:none;border-radius:4px;cursor:pointer;font-weight:500;';
      selectButton.addEventListener('click', function() {
        selectMolecule(molecule.smiles, card);
      });
      
      // Add delete button
      const deleteButton = document.createElement('button');
      deleteButton.textContent = 'DELETE';
      deleteButton.style.cssText = 'flex:1;padding:8px;background:#f44336;color:white;border:none;border-radius:4px;cursor:pointer;font-weight:500;';
      deleteButton.addEventListener('click', function() {
        deleteMolecule(molecule.id);
      });
      
      buttonsContainer.appendChild(selectButton);
      buttonsContainer.appendChild(deleteButton);
      card.appendChild(buttonsContainer);
      
      container.appendChild(card);
      
      // Add visualization
      const previewContainer = card.querySelector('.molecule-preview');
      if (previewContainer && typeof createMoleculeVisualization === 'function') {
        createMoleculeVisualization(previewContainer, molecule.smiles);
      }
    });
  } catch (error) {
    console.error('Error updating saved molecules tab:', error);
  }
}

// Delete a saved molecule
function deleteMolecule(id) {
  try {
    // Get saved molecules
    const savedMolecules = JSON.parse(localStorage.getItem('savedMolecules') || '[]');
    
    // Find and remove the molecule
    const index = savedMolecules.findIndex(m => m.id === id);
    if (index !== -1) {
      savedMolecules.splice(index, 1);
      localStorage.setItem('savedMolecules', JSON.stringify(savedMolecules));
      
      // Update UI
      updateSavedMoleculesTab();
      
      // Show feedback
      showToast('Molecule deleted');
    }
  } catch (error) {
    console.error('Error deleting molecule:', error);
    showToast('Error deleting molecule', 'error');
  }
}

// Reuse the isSmilesFormat and extractSmiles functions from other scripts
function isSmilesFormat(text) {
  if (!text) return false;
  text = text.trim();
  
  // If labeled as SMILES
  if (text.includes('SMILES:')) {
    const match = text.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
    return match ? true : false;
  }
  
  // Check if it follows SMILES format (basic validation)
  if (text.length > 5 && text.length < 200 && 
      /^[A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+$/.test(text) &&
      /[C|c|N|n|O|o]/.test(text)) {
    return true;
  }
  
  return false;
}

function extractSmiles(text) {
  if (!text) return '';
  text = text.trim();
  
  // If labeled as SMILES
  if (text.includes('SMILES:')) {
    const match = text.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
    return match ? match[1].trim() : text;
  }
  
  return text;
}

// Toast notification function (shared across scripts)
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

// Function to replace AI with Claude across the interface
function replaceAIWithClaude() {
  // Replace text content in elements
  const elements = document.querySelectorAll('h1, h2, h3, h4, h5, h6, p, span, div, button, a, label');
  
  elements.forEach(element => {
    // Skip elements with children (to avoid duplicating replacements)
    if (element.children.length > 0) return;
    
    // Replace text content
    if (element.textContent.includes('AI')) {
      element.textContent = element.textContent.replace(/\bAI\b/g, 'Claude');
    }
  });
  
  // Replace class names and IDs
  document.querySelectorAll('[class*="ai"], [id*="ai"]').forEach(element => {
    // Get all classes
    const classes = Array.from(element.classList);
    classes.forEach(className => {
      if (className.includes('ai')) {
        const newClass = className.replace(/ai/gi, 'claude');
        element.classList.add(newClass);
      }
    });
    
    // Update ID if needed
    if (element.id && element.id.includes('ai')) {
      element.id = element.id.replace(/ai/gi, 'claude');
    }
  });
  
  // Update avatar letters from AI to C
  document.querySelectorAll('.MuiAvatar-root, [class*="avatar"]').forEach(avatar => {
    if (avatar.textContent === 'AI') {
      avatar.textContent = 'C';
    }
  });
}
