// Enhanced Fix for 3D molecule visualization

// Wait for document to be ready
document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying aggressive molecule visualization fixes...');
  
  // Apply fixes immediately
  fixMoleculeVisualization();
  
  // First 10 seconds - check very frequently (every 200ms)
  for (let i = 1; i <= 50; i++) {
    setTimeout(() => fixMoleculeVisualization(), i * 200);
  }
  
  // Then check periodically
  setInterval(fixMoleculeVisualization, 1000);
  
  // Set up an observer to detect when new molecules are added
  setupMoleculeObserver();
});

// Set up mutation observer to watch for new molecules
function setupMoleculeObserver() {
  console.log('Setting up molecule observer');
  
  // Create an observer that watches for changes in the DOM
  const observer = new MutationObserver(function(mutations) {
    let shouldCheckVisualization = false;
    
    // Check if any of the mutations involve adding nodes
    mutations.forEach(function(mutation) {
      if (mutation.addedNodes.length > 0) {
        // Check if any added nodes contain molecule-related content
        Array.from(mutation.addedNodes).forEach(function(node) {
          if (node.nodeType === 1) { // Element node
            const innerHTML = node.innerHTML || '';
            if (innerHTML.includes('SMILES') || 
                innerHTML.includes('molecule') || 
                innerHTML.includes('Molecule') || 
                innerHTML.includes('structure')) {
              shouldCheckVisualization = true;
            }
          }
        });
      }
    });
    
    if (shouldCheckVisualization) {
      console.log('Detected new molecule content, applying visualization fixes');
      fixMoleculeVisualization();
      
      // Check again after short delay to ensure everything is rendered
      setTimeout(fixMoleculeVisualization, 500);
      setTimeout(fixMoleculeVisualization, 1000);
    }
  });
  
  // Start observing the document with the configured parameters
  observer.observe(document.body, { childList: true, subtree: true });
}

// Main visualization fix function
function fixMoleculeVisualization() {
  // First check if success message is shown but no molecules are visible
  const successMessage = document.querySelector('.success-message, [class*="success"]');
  if (successMessage) {
    console.log('Success message found - ensuring molecules are displayed');
    showAllMoleculeContent();
  }
  
  // Find all possible molecule containers that might need visualization
  const viewerContainers = document.querySelectorAll(
    '.viewerContainer, [class*="moleculeViewer"], [class*="Viewer3D"], .molecule-preview, ' + 
    '.model-container, .visualization, [class*="model"], [class*="3d"], [class*="viewer"]'
  );
  
  console.log('Found ' + viewerContainers.length + ' potential molecule containers');
  
  // Process each container
  viewerContainers.forEach(function(container) {
    // Skip if already successfully processed and has visible content
    if (container.dataset.processed === 'true' && container.innerHTML.length > 100) return;
    
    // Find SMILES string
    let smiles = findSmilesNearContainer(container);
    if (!smiles) {
      // If no SMILES found directly, search more aggressively
      smiles = findSmilesAnywhere();
      if (!smiles) return;
    }
    
    console.log('Found SMILES to visualize:', smiles.substring(0, 20) + '...');
    
    // Create a visualization
    createMoleculeVisualization(container, smiles);
    
    // Mark as processed
    container.dataset.processed = 'true';
  });

  // Fix empty "Generated Molecules" section
  fixGeneratedMolecules();
  
  // Fix any molecule cards that are missing info
  enhanceMoleculeCards();
}

// Find SMILES string near a container
function findSmilesNearContainer(container) {
  // Look in the container itself
  let smilesElement = container.querySelector('.smiles, [class*="smiles"], pre, code');
  if (smilesElement && smilesElement.textContent && isSmilesFormat(smilesElement.textContent)) {
    return smilesElement.textContent.trim();
  }
  
  // Look in parent container
  const parent = container.closest('.card, [class*="Card"], [class*="molecule"]') || container.parentElement;
  if (!parent) return null;
  
  // Check for SMILES in various elements
  const possibleElements = parent.querySelectorAll('pre, code, .smiles, [class*="smiles"], p, div');
  for (const element of possibleElements) {
    if (element.textContent && isSmilesFormat(element.textContent)) {
      return element.textContent.trim();
    }
  }
  
  return null;
}

// Check if a string looks like SMILES format
function isSmilesFormat(text) {
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

// Extract SMILES from text
function extractSmiles(text) {
  text = text.trim();
  
  // If labeled as SMILES
  if (text.includes('SMILES:')) {
    const match = text.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
    return match ? match[1].trim() : text;
  }
  
  return text;
}

// Create a 3D visualization 
function createMoleculeVisualization(container, smilesText) {
  try {
    // Clean up SMILES format
    const smiles = extractSmiles(smilesText);
    
    // Create a unique ID
    const viewerId = 'mol-viewer-' + Math.random().toString(36).substring(2, 9);
    
    // Create a container for the viewer
    container.innerHTML = `
      <div id="${viewerId}" style="width:100%;height:100%;min-height:200px;position:relative;">
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
    
    // Load 3Dmol.js if needed
    if (typeof $3Dmol === 'undefined') {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.org/build/3Dmol-min.js';
      document.head.appendChild(script);
      
      script.onload = function() {
        create3DMolVisualization(smiles, viewerElement);
      };
    } else {
      create3DMolVisualization(smiles, viewerElement);
    }
  } catch (error) {
    console.error('Error creating visualization:', error);
    container.innerHTML = `<div style="color:red;padding:10px;text-align:center;">Error creating molecule visualization</div>`;
  }
}

// Create the actual 3D visualization
function create3DMolVisualization(smiles, viewerElement) {
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
    viewerElement.innerHTML = `<div style="color:red;padding:10px;text-align:center;">Error: ${error.message}</div>`;
  });
}

// Find SMILES strings anywhere in the document
function findSmilesAnywhere() {
  console.log('Looking for SMILES strings anywhere in the document');
  
  // Look in thinking process sections first
  const thinkingProcess = document.querySelector('.ThinkingProcess, [class*="thinkingProcess"], [class*="thinking"]');
  if (thinkingProcess) {
    const preElements = thinkingProcess.querySelectorAll('pre, code');
    for (const element of preElements) {
      if (element.textContent && isSmilesFormat(element.textContent)) {
        return element.textContent.trim();
      }
    }
    
    // Look for text that contains "SMILES" label
    const allElements = thinkingProcess.querySelectorAll('*');
    for (const element of allElements) {
      if (element.textContent && element.textContent.includes('SMILES:')) {
        const match = element.textContent.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
        if (match) return match[1].trim();
      }
    }
  }
  
  // Look in any element that contains "SMILES"
  const allElements = document.querySelectorAll('*');
  for (const element of allElements) {
    if (element.textContent && element.textContent.includes('SMILES:')) {
      const match = element.textContent.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
      if (match) return match[1].trim();
    }
  }
  
  // As a last resort, look for anything that looks like a SMILES string
  const textElements = document.querySelectorAll('p, pre, code, span, div');
  for (const element of textElements) {
    if (element.textContent && isSmilesFormat(element.textContent)) {
      return element.textContent.trim();
    }
  }
  
  return null;
}

// Show all molecule content in the document
function showAllMoleculeContent() {
  console.log('Forcibly showing all molecule content');
  
  // Find all molecule-related elements
  const selectors = [
    '.molecule-card', '.molecule', '[class*="molecule"]',
    '.smiles-string', '.smiles', '[class*="smiles"]',
    '.molecule-details', '.properties', '.molecular-properties',
    '.structure', '.model', '.visualization',
    '.generated-molecules', '[class*="generated"]'
  ];
  
  selectors.forEach(selector => {
    const elements = document.querySelectorAll(selector);
    console.log(`Found ${elements.length} elements matching ${selector}`);
    
    elements.forEach(element => {
      // Force display
      element.style.display = 'block';
      element.style.visibility = 'visible';
      element.style.opacity = '1';
      
      // If parent is hidden, show it too
      let parent = element.parentElement;
      for (let i = 0; i < 5 && parent; i++) { // Check up to 5 levels up
        parent.style.display = 'block';
        parent.style.visibility = 'visible';
        parent = parent.parentElement;
      }
    });
  });
  
  // Also find elements containing keyword-based content and show them
  const allElements = document.querySelectorAll('div, section, article');
  allElements.forEach(element => {
    if (element.innerHTML && (
        element.innerHTML.includes('SMILES') ||
        element.innerHTML.includes('molecule') ||
        element.innerHTML.includes('Molecule') ||
        element.innerHTML.includes('structure')
    )) {
      element.style.display = 'block';
      element.style.visibility = 'visible';
      element.style.opacity = '1';
    }
  });
  
  // Add CSS to ensure visibility
  if (!document.querySelector('#molecule-visibility-fix')) {
    const styleEl = document.createElement('style');
    styleEl.id = 'molecule-visibility-fix';
    styleEl.textContent = `
      /* Force molecule visibility */
      .molecule-card, .molecule, [class*="molecule"],
      .smiles-string, .smiles, [class*="smiles"],
      .molecule-details, .properties, .molecular-properties,
      .structure, .model, .visualization,
      .generated-molecules, [class*="generated"] {
        display: block !important;
        visibility: visible !important;
        opacity: 1 !important;
      }
    `;
    document.head.appendChild(styleEl);
  }
}

// Enhance molecule cards with better styling and data
function enhanceMoleculeCards() {
  // Find all molecule cards
  const moleculeCards = document.querySelectorAll('.molecule-card, [class*="molecule"][class*="card"]');
  console.log(`Enhancing ${moleculeCards.length} molecule cards`);
  
  moleculeCards.forEach(card => {
    // Skip if already enhanced
    if (card.dataset.enhanced === 'true') return;
    
    // Find SMILES string
    let smilesElement = card.querySelector('.smiles-string, [class*="smiles"]');
    let smiles = smilesElement ? smilesElement.textContent.trim() : null;
    
    if (!smiles) {
      // Look for SMILES in the card content
      const content = card.innerHTML;
      if (content.includes('SMILES:')) {
        const match = content.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
        if (match) smiles = match[1].trim();
      }
    }
    
    if (smiles) {
      // Add styling to the card if not already styled
      if (!card.style.border) {
        card.style.border = '1px solid rgba(0,0,0,0.12)';
        card.style.borderRadius = '4px';
        card.style.padding = '16px';
        card.style.margin = '8px 0';
        card.style.backgroundColor = 'white';
        card.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
      }
      
      // Make sure the card has a visualization container
      let visualizationContainer = card.querySelector('.model-container, .visualization, [class*="model"]');
      if (!visualizationContainer) {
        visualizationContainer = document.createElement('div');
        visualizationContainer.className = 'model-container';
        visualizationContainer.style.height = '300px';
        visualizationContainer.style.width = '100%';
        visualizationContainer.style.marginTop = '16px';
        visualizationContainer.style.marginBottom = '16px';
        visualizationContainer.style.backgroundColor = '#f5f5f5';
        visualizationContainer.style.borderRadius = '4px';
        
        // Find a good insertion point
        const smilesContainer = card.querySelector('.smiles-string, [class*="smiles"]');
        if (smilesContainer) {
          smilesContainer.parentNode.insertBefore(visualizationContainer, smilesContainer.nextSibling);
        } else {
          card.appendChild(visualizationContainer);
        }
        
        // Create the visualization
        createMoleculeVisualization(visualizationContainer, smiles);
      }
      
      // Make sure SMILES string is displayed properly
      if (smilesElement) {
        smilesElement.style.fontFamily = 'monospace';
        smilesElement.style.wordBreak = 'break-all';
        smilesElement.style.backgroundColor = '#f5f5f5';
        smilesElement.style.padding = '8px';
        smilesElement.style.borderRadius = '4px';
      }
    }
    
    // Mark as enhanced
    card.dataset.enhanced = 'true';
  });
}

// Fix Generated Molecules section
function fixGeneratedMolecules() {
  const generatedMoleculesSection = document.querySelector('.generated-molecules, .GeneratedMolecules, [class*="generatedMolecules"]');
  if (!generatedMoleculesSection) {
    // Try to find a container that might be the generated molecules section
    const possibleSections = document.querySelectorAll('h2, h3, h4');
    for (const heading of possibleSections) {
      if (heading.textContent.includes('Generated Molecules') || 
          heading.textContent.includes('Molecule') || 
          heading.textContent.includes('Results')) {
        const parent = heading.parentElement;
        if (parent) {
          console.log('Found potential generated molecules container:', heading.textContent);
          parent.style.display = 'block';
          parent.style.visibility = 'visible';
          
          // Check if it's empty and needs populating
          const hasContent = parent.querySelectorAll('.molecule-card, [class*="molecule"]').length > 0;
          if (!hasContent) {
            populateWithMolecules(parent);
          }
        }
      }
    }
    return;
  }
  
  console.log('Found generated molecules section');
  generatedMoleculesSection.style.display = 'block';
  generatedMoleculesSection.style.visibility = 'visible';
  
  // Check if empty
  const emptyMessage = generatedMoleculesSection.querySelector('p, div');
  const hasCards = generatedMoleculesSection.querySelectorAll('.molecule-card, [class*="molecule"]').length > 0;
  
  if ((!hasCards && emptyMessage && emptyMessage.textContent.includes('No')) || 
      generatedMoleculesSection.innerHTML.trim() === '') {
    console.log('Empty generated molecules section, attempting to populate');
    populateWithMolecules(generatedMoleculesSection);
  }
}

// Populate a container with molecule cards from SMILES strings found elsewhere
function populateWithMolecules(container) {
  // Find SMILES strings in ThinkingProcess
  const thinkingProcess = document.querySelector('.ThinkingProcess, [class*="thinking"], [class*="response"]');
  if (!thinkingProcess) {
    console.log('No thinking process found to extract molecules from');
    return;
  }
  
  // Try to find SMILES strings
  const smilesElements = Array.from(thinkingProcess.querySelectorAll('pre, code, [class*="smiles"]'))
    .filter(element => isSmilesFormat(element.textContent));
  
  console.log(`Found ${smilesElements.length} SMILES strings in thinking process`);
  
  if (smilesElements.length === 0) {
    // Look for text blocks that might contain SMILES
    const textBlocks = thinkingProcess.querySelectorAll('p, div, span');
    for (const block of textBlocks) {
      if (block.textContent && block.textContent.includes('SMILES:')) {
        const match = block.textContent.match(/SMILES:?\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
        if (match) {
          smilesElements.push({
            textContent: match[1].trim()
          });
        }
      }
    }
  }
  
  if (smilesElements.length === 0) {
    console.log('No SMILES strings found');
    return;
  }
  
  // Clear empty message if present
  if (container.querySelector('p')?.textContent.includes('No')) {
    container.innerHTML = '';
  }
  
  // Create molecule cards
  const moleculesContainer = document.createElement('div');
  moleculesContainer.className = 'molecules-grid';
  moleculesContainer.style.display = 'grid';
  moleculesContainer.style.gridTemplateColumns = 'repeat(auto-fill, minmax(300px, 1fr))';
  moleculesContainer.style.gap = '16px';
  moleculesContainer.style.marginTop = '16px';
  
  smilesElements.forEach((element, index) => {
    const smiles = extractSmiles(element.textContent);
      
    const card = document.createElement('div');
    card.className = 'molecule-card';
    card.style.width = '300px';
    card.style.border = '1px solid #e0e0e0';
    card.style.borderRadius = '4px';
    card.style.overflow = 'hidden';
    
    card.innerHTML = `
      <div style="padding:12px;background-color:#f5f5f5;">
        <h3 style="margin:0;">Molecule ${index + 1}</h3>
      </div>
      <div class="molecule-preview" style="height:200px;"></div>
      <div style="padding:12px;">
        <pre style="overflow:auto;background:#f9f9f9;padding:8px;font-size:12px;">${smiles}</pre>
        <div style="display:flex;gap:8px;margin-top:12px;">
          <button class="select-btn" style="flex:1;padding:8px;background:#2196f3;color:white;border:none;border-radius:4px;cursor:pointer;">SELECT</button>
          <button class="save-btn" style="flex:1;padding:8px;background:#4caf50;color:white;border:none;border-radius:4px;cursor:pointer;">SAVE</button>
        </div>
      </div>
    `;
    
    moleculesContainer.appendChild(card);
    
    // Find molecule preview element and visualize molecule
    const previewEl = card.querySelector('.molecule-preview');
    if (previewEl) {
      createMoleculeVisualization(previewEl, smiles);
    }
  });
  
  // Replace empty message with molecules
  if (container.innerHTML.includes('No') || container.innerHTML.trim() === '') {
    container.innerHTML = '';
  }
  container.appendChild(moleculesContainer);
}
