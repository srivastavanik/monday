// Enhanced molecule visualization and display restoration script

document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying molecule visualization restoration...');
  
  // Apply fixes immediately and repeatedly with increasing frequency
  restoreAllContent();
  
  // Initial checks every 100ms for first 5 seconds
  let counter = 0;
  const fastInterval = setInterval(() => {
    restoreAllContent();
    counter++;
    if (counter >= 50) { // 5 seconds
      clearInterval(fastInterval);
    }
  }, 100);
  
  // Then check every second indefinitely
  setInterval(restoreAllContent, 1000);
});

function restoreAllContent() {
  try {
    // Process success message first - often indicates molecules were generated
    const successMsg = document.querySelector('.success-message, [class*="success"], .success');
    if (successMsg) {
      console.log('Success message found, ensuring molecules are displayed');
      findAndShowMolecules();
    }
    
    // 1. Show generated molecules container with high specificity
    const generatedContainers = document.querySelectorAll('.generated-molecules, [class*="generated"], [class*="candidates"]');
    generatedContainers.forEach(container => {
      if (container) {
        container.style.display = 'block';
        container.style.visibility = 'visible';
        container.style.maxHeight = 'none';
        container.style.overflow = 'visible';
        container.style.opacity = '1';
      }
    });
    
    // 2. Find and show all molecule-related elements
    showMoleculeElements();
    
    // 3. Chat restoration
    const chatElements = document.querySelectorAll('.message-container, .chat-container, .claude-response, .ai-message');
    chatElements.forEach(element => {
      if (element) {
        element.style.display = 'block';
        element.style.visibility = 'visible';
        element.style.maxHeight = 'none';
        element.style.overflow = 'visible';
      }
    });
    
    // 4. Add comprehensive CSS to ensure correct display
    ensureDisplayStyles();
    
    // 5. Attempt to render any 3D models that haven't been rendered
    tryRenderMoleculeModels();
    
  } catch (err) {
    console.error('Error in display restoration:', err);
  }
}

function showMoleculeElements() {
  // Target all possible molecule-related elements with high specificity
  const selectors = [
    '.molecule-card', '[class*="molecule"]', '[class*="structure"]',
    '.smiles-string', '.molecule-details', '.molecule-name',
    '.molecule-actions', '.save-molecule', '.molecule-properties',
    '.model-container', '.visualization', '.3d-model',
    '.property-table', '.molecular-properties', '.rdkit-canvas'
  ];
  
  selectors.forEach(selector => {
    const elements = document.querySelectorAll(selector);
    elements.forEach(element => {
      if (element) {
        // Force display
        element.style.display = 'block';
        element.style.visibility = 'visible';
        element.style.opacity = '1';
        
        // If it's a container, ensure content is visible
        if (element.children && element.children.length > 0) {
          Array.from(element.children).forEach(child => {
            child.style.display = 'block';
            child.style.visibility = 'visible';
          });
        }
      }
    });
  });
}

function findAndShowMolecules() {
  // Find any divs that might contain molecule information
  const allDivs = document.querySelectorAll('div');
  allDivs.forEach(div => {
    // Check if div has molecule content
    if (div.innerHTML && (
        div.innerHTML.includes('SMILES') || 
        div.innerHTML.includes('Molecule') || 
        div.innerHTML.includes('molecule') ||
        div.innerHTML.includes('structure') ||
        div.innerHTML.includes('molecular') ||
        div.innerHTML.includes('3D model') ||
        div.innerHTML.includes('visualize') ||
        div.innerHTML.toLowerCase().includes('properties'))) {
      
      // Show this container and all its children
      div.style.display = 'block';
      div.style.visibility = 'visible';
      div.style.opacity = '1';
      div.style.maxHeight = 'none';
      
      // Also show all children
      if (div.children && div.children.length > 0) {
        Array.from(div.children).forEach(child => {
          child.style.display = 'block';
          child.style.visibility = 'visible';
        });
      }
    }
  });
}

function tryRenderMoleculeModels() {
  // Look for molecule visualization containers that might need rendering
  const modelContainers = document.querySelectorAll('.model-container, .visualization, .3d-model, [id*="molecule"], [class*="3d"]');
  
  modelContainers.forEach(container => {
    // If container is empty or has minimal content, try to trigger rendering
    if (container && (!container.innerHTML || container.innerHTML.trim().length < 50)) {
      // Try to find SMILES data
      const parentCard = findParentByClass(container, 'molecule-card') || 
                         findParentByClass(container, 'molecule') ||
                         container.closest('[class*="molecule"]');
      
      if (parentCard) {
        // Find SMILES string in this card
        const smilesElement = parentCard.querySelector('.smiles-string') || 
                             parentCard.querySelector('[class*="smiles"]');
        
        if (smilesElement && smilesElement.textContent) {
          const smiles = smilesElement.textContent.trim();
          if (smiles && smiles.length > 0) {
            console.log('Found empty model container with SMILES, triggering render:', smiles);
            // Dispatch custom event to trigger rendering
            const event = new CustomEvent('render-molecule', { 
              detail: { 
                container: container,
                smiles: smiles
              }
            });
            document.dispatchEvent(event);
          }
        }
      }
    }
  });
}

function findParentByClass(element, className) {
  let current = element;
  while (current) {
    if (current.classList && current.classList.contains(className)) {
      return current;
    }
    current = current.parentElement;
  }
  return null;
}

function ensureDisplayStyles() {
  if (!document.querySelector('#restore-styles')) {
    const styleEl = document.createElement('style');
    styleEl.id = 'restore-styles';
    styleEl.innerHTML = `
      /* Comprehensive molecule display fix */
      .molecule-card, [class*="molecule"], [class*="structure"],
      .smiles-string, .molecule-details, .molecule-name,
      .molecule-actions, .save-molecule, .molecule-properties,
      .model-container, .visualization, .3d-model,
      .property-table, .molecular-properties, .rdkit-canvas,
      .generated-molecules, [class*="generated"], [class*="candidates"] {
        display: block !important;
        visibility: visible !important;
        opacity: 1 !important;
        max-height: none !important;
      }
      
      /* Make molecule cards more prominent */
      .molecule-card {
        border: 1px solid rgba(0, 0, 0, 0.12) !important;
        border-radius: 4px !important;
        margin-bottom: 20px !important;
        padding: 16px !important;
        background-color: white !important;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1) !important;
      }
      
      /* Make model containers at least 300px height */
      .model-container, .visualization, .3d-model {
        min-height: 300px !important;
      }
      
      /* Style SMILES strings */
      .smiles-string {
        font-family: monospace !important;
        word-break: break-all !important;
        background-color: #f5f5f5 !important;
        padding: 8px !important;
        border-radius: 4px !important;
        margin: 8px 0 !important;
      }
      
      /* Ensure chat containers are visible */
      .message-container, .chat-container, .claude-response, .ai-message {
        display: block !important;
        visibility: visible !important;
        max-height: none !important;
        overflow: visible !important;
      }
      
      /* Fix the height of the status message */
      .claude-is-designing {
        max-height: 100px !important;
        overflow: hidden !important;
      }
    `;
    document.head.appendChild(styleEl);
  }
}

// Listen for custom render events
document.addEventListener('render-molecule', function(e) {
  if (e.detail && e.detail.container && e.detail.smiles) {
    console.log('Received render event for:', e.detail.smiles);
    try {
      // If 3Dmol is available, try to render
      if (window.$3Dmol) {
        const container = e.detail.container;
        const smiles = e.detail.smiles;
        
        // Set a unique ID if not present
        if (!container.id) {
          container.id = 'molecule-viz-' + Math.random().toString(36).substring(2, 9);
        }
        
        // Ensure container has height
        container.style.height = '300px';
        container.style.width = '100%';
        container.style.position = 'relative';
        
        // Create viewer
        const viewer = window.$3Dmol.createViewer(container.id, {
          backgroundColor: 'white'
        });
        
        if (viewer) {
          // Try to add the molecule from SMILES
          const smilesParser = new window.$3Dmol.SmilesParser();
          smilesParser.addMolData(viewer, smiles, 'smi', {}, function() {
            viewer.setStyle({}, {stick: {}}); 
            viewer.zoomTo();
            viewer.render();
            console.log('Rendered 3D model for', smiles);
          });
        }
      }
    } catch (err) {
      console.error('Error rendering molecule:', err);
    }
  }
});
