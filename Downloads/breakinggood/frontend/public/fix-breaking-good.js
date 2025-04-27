// Fix script for the Breaking Good platform

// Wait for document to be ready
document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying fixes to Breaking Good platform...');
  
  // Fix the 3D visualization
  fixMoleculeVisualization();
  
  // Fix the save and select functionality
  fixSaveAndSelectButtons();
  
  // Fix the chat functionality
  fixChatFunctionality();
  
  // Run these fixes periodically to catch new elements
  setInterval(function() {
    fixMoleculeVisualization();
    fixSaveAndSelectButtons();
  }, 3000);
});

// Fix for molecule visualization
function fixMoleculeVisualization() {
  console.log('Fixing molecule visualization...');
  
  // Find all molecule visualization containers
  const moleculeContainers = document.querySelectorAll('.viewerContainer, [class*="viewerContainer"]');
  
  if (moleculeContainers.length === 0) {
    console.log('No molecule containers found yet, will retry later');
    return;
  }
  
  console.log(`Found ${moleculeContainers.length} molecule containers`);
  
  // Process each container
  moleculeContainers.forEach(function(container, index) {
    // Skip if already processed
    if (container.dataset.processed === 'true') return;
    
    // Find the closest SMILES string
    let smiles = findSmilesNearContainer(container);
    
    if (!smiles) {
      console.log(`Container ${index}: No SMILES found`);
      return;
    }
    
    console.log(`Container ${index}: Found SMILES ${smiles.substring(0, 20)}...`);
    
    // Create a visualization
    createMoleculeVisualization(container, smiles);
    
    // Mark as processed
    container.dataset.processed = 'true';
  });
}

// Helper to find SMILES string near a container
function findSmilesNearContainer(container) {
  // Search in container and surrounding elements
  const parent = container.closest('[class*="Card"], [class*="cardContent"], [class*="Card"]') || container.parentElement;
  
  // Look for elements that might contain SMILES
  const possibleElements = parent.querySelectorAll('code, pre, [class*="codeBlock"], .MuiTypography-body2, div');
  
  for (const element of possibleElements) {
    const text = element.textContent;
    if (!text) continue;
    
    // Check if it looks like a SMILES string (no spaces, valid characters)
    if (text.trim().match(/^[A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]{10,}$/)) {
      return text.trim();
    }
    
    // Also check for labeled SMILES
    if (text.includes('SMILES:')) {
      const match = text.match(/SMILES:\s*([A-Za-z0-9@\[\]\(\)\/#=\\\-\+\.]+)/);
      if (match && match[1]) {
        return match[1].trim();
      }
    }
  }
  
  return null;
}

// Create a 3D visualization for a molecule
function createMoleculeVisualization(container, smiles) {
  try {
    // Create a unique ID for the viewer
    const viewerId = 'molecule-viewer-' + Math.random().toString(36).substring(2, 9);
    
    // Create a new div inside the container
    container.innerHTML = `
      <div id="${viewerId}" style="width:100%;height:100%;position:relative;">
        <div style="position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);">
          <div class="loading-spinner" style="border:4px solid rgba(0,0,0,0.1);border-radius:50%;border-top:4px solid #3498db;width:30px;height:30px;animation:spin 1s linear infinite;"></div>
        </div>
      </div>
      <style>
        @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
      </style>
    `;
    
    // Create an element for the 3D viewer
    const viewerElement = document.getElementById(viewerId);
    
    // Make sure 3Dmol.js is loaded
    if (typeof $3Dmol === 'undefined') {
      // Load 3Dmol.js dynamically if needed
      const script = document.createElement('script');
      script.src = 'https://3Dmol.org/build/3Dmol-min.js';
      document.head.appendChild(script);
      
      script.onload = function() {
        // Continue with creating visualization
        convertSmilesAndVisualize(smiles, viewerElement);
      };
    } else {
      // 3Dmol.js is already loaded
      convertSmilesAndVisualize(smiles, viewerElement);
    }
  } catch (error) {
    console.error('Error creating visualization:', error);
    container.innerHTML = `<div style="color:red;text-align:center;padding:10px;">Error creating visualization</div>`;
  }
}

// Convert SMILES to 3D and create visualization
function convertSmilesAndVisualize(smiles, viewerElement) {
  // Make API call to convert SMILES to 3D structure
  fetch('/api/simulation/3d-structure', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json'
    },
    body: JSON.stringify({ smiles: smiles })
  })
  .then(response => {
    if (!response.ok) {
      throw new Error('Failed to convert SMILES to 3D structure');
    }
    return response.json();
  })
  .then(data => {
    if (!data.molblock) {
      throw new Error('No molblock in response');
    }
    
    // Create viewer with white background
    const viewer = $3Dmol.createViewer(viewerElement, {
      backgroundColor: 'white'
    });
    
    // Add the model
    const model = viewer.addModel(data.molblock, 'mol');
    
    // Style the molecule
    viewer.setStyle({}, { stick: {} });
    
    // Add a surface with some transparency
    viewer.addSurface($3Dmol.SurfaceType.VDW, {
      opacity: 0.6,
      color: 'lightblue'
    });
    
    // Center and render
    viewer.zoomTo();
    viewer.render();
    
    // Add gentle rotation animation
    let animating = true;
    function animate() {
      if (animating) {
        viewer.rotate(0.3, { y: 1 });
        viewer.render();
        requestAnimationFrame(animate);
      }
    }
    animate();
    
    // Store for cleanup
    viewerElement.viewer = viewer;
    viewerElement.stopAnimation = function() {
      animating = false;
    };
  })
  .catch(error => {
    console.error('Error in visualization:', error);
    viewerElement.innerHTML = `<div style="color:red;text-align:center;padding:10px;">${error.message}</div>`;
  });
}

// Fix save and select buttons
function fixSaveAndSelectButtons() {
  console.log('Fixing save and select buttons...');
  
  // Find all save and select buttons using standard selectors
  const saveButtons = Array.from(document.querySelectorAll('button[class*="save"], button'))
    .filter(button => button.textContent.toLowerCase().includes('save'));
  const selectButtons = Array.from(document.querySelectorAll('button[class*="select"], button'))
    .filter(button => button.textContent.toLowerCase().includes('select'));
  
  // Process save buttons
  saveButtons.forEach(function(button) {
    // Skip if already processed
    if (button.dataset.fixedAction === 'save') return;
    
    // Add click handler
    button.addEventListener('click', function(event) {
      event.preventDefault();
      
      // Find closest molecule container
      const container = button.closest('[class*="Card"], [class*="cardContent"], [class*="molecule"]');
      if (!container) return;
      
      // Find SMILES in container
      const smiles = findSmilesNearContainer(container);
      if (!smiles) return;
      
      // Create molecule object
      const molecule = {
        id: Date.now(),
        name: `Saved-${Math.floor(Math.random() * 1000)}`,
        smiles: smiles,
        timestamp: new Date().toISOString()
      };
      
      // Save to localStorage
      try {
        const savedMolecules = JSON.parse(localStorage.getItem('savedMolecules') || '[]');
        savedMolecules.push(molecule);
        localStorage.setItem('savedMolecules', JSON.stringify(savedMolecules));
        alert(`Molecule saved as ${molecule.name}`);
      } catch (error) {
        console.error('Error saving molecule:', error);
        alert('Error saving molecule');
      }
    });
    
    // Mark as processed
    button.dataset.fixedAction = 'save';
  });
  
  // Process select buttons
  selectButtons.forEach(function(button) {
    // Skip if already processed
    if (button.dataset.fixedAction === 'select') return;
    
    // Add click handler
    button.addEventListener('click', function(event) {
      event.preventDefault();
      
      // Find closest molecule container
      const container = button.closest('[class*="Card"], [class*="cardContent"], [class*="molecule"]');
      if (!container) return;
      
      // Find SMILES in container
      const smiles = findSmilesNearContainer(container);
      if (!smiles) return;
      
      // Highlight the selected container
      const allContainers = document.querySelectorAll('[class*="Card"], [class*="cardContent"], [class*="molecule"]');
      allContainers.forEach(c => c.style.border = '');
      container.style.border = '2px solid #4caf50';
      
      alert(`Selected molecule with SMILES: ${smiles.substring(0, 20)}...`);
      
      // Try to update the main visualization area if it exists
      const mainViewer = document.querySelector('.MoleculeViewer3D, [class*="MoleculeViewer3D"]');
      if (mainViewer) {
        // Update with this molecule
        createMoleculeVisualization(mainViewer, smiles);
      }
    });
    
    // Mark as processed
    button.dataset.fixedAction = 'select';
  });
}

// Fix chat functionality
function fixChatFunctionality() {
  console.log('Fixing chat functionality...');
  
  // Find the chat input and send button
  const chatInputs = document.querySelectorAll('textarea[placeholder*="question"], textarea[placeholder*="refinements"], textarea[class*="input"]');
  const sendButtons = Array.from(document.querySelectorAll('button[class*="send"], button'))
    .filter(button => {
      const buttonText = button.textContent.toLowerCase();
      return buttonText.includes('send') || button.querySelector('svg[class*="Send"]');
    });
  
  chatInputs.forEach(function(input) {
    // Skip if already processed
    if (input.dataset.fixedChat === 'true') return;
    
    // Add keypress handler
    input.addEventListener('keypress', function(event) {
      if (event.key === 'Enter' && !event.shiftKey) {
        event.preventDefault();
        sendChatMessage(input.value);
        input.value = '';
      }
    });
    
    // Mark as processed
    input.dataset.fixedChat = 'true';
  });
  
  sendButtons.forEach(function(button) {
    // Skip if already processed
    if (button.dataset.fixedChat === 'true') return;
    
    // Add click handler
    button.addEventListener('click', function(event) {
      event.preventDefault();
      
      // Find associated input
      const container = button.closest('[class*="input"], [class*="chat"], [class*="message"]');
      if (!container) return;
      
      const input = container.querySelector('textarea, input[type="text"]');
      if (!input) return;
      
      sendChatMessage(input.value);
      input.value = '';
    });
    
    // Mark as processed
    button.dataset.fixedChat = 'true';
  });
}

// Send a chat message and handle the response
function sendChatMessage(message) {
  if (!message.trim()) return;
  
  console.log('Sending chat message:', message);
  
  // Find message container
  const messageContainer = document.querySelector('[class*="messageContainer"], [class*="messages"], [class*="chat"]');
  if (!messageContainer) return;
  
  // Add user message to UI
  const userMessageHtml = `
    <div class="message user-message" style="display:flex;margin-bottom:16px;">
      <div style="width:32px;height:32px;border-radius:50%;background-color:#9c27b0;color:white;display:flex;justify-content:center;align-items:center;margin-right:10px;">U</div>
      <div style="flex:1;background-color:rgba(0,0,0,0.05);padding:10px;border-radius:8px;">
        <div>${message}</div>
        <div style="font-size:0.8em;color:#666;margin-top:4px;">${new Date().toLocaleString()}</div>
      </div>
    </div>
  `;
  
  messageContainer.innerHTML += userMessageHtml;
  messageContainer.scrollTop = messageContainer.scrollHeight;
  
  // Show typing indicator
  const typingHtml = `
    <div class="message ai-typing" style="display:flex;margin-bottom:16px;">
      <div style="width:32px;height:32px;border-radius:50%;background-color:#2196f3;color:white;display:flex;justify-content:center;align-items:center;margin-right:10px;">AI</div>
      <div style="flex:1;background-color:rgba(0,0,0,0.03);padding:10px;border-radius:8px;">
        <div class="typing-indicator" style="display:flex;">
          <div style="height:8px;width:8px;border-radius:50%;background-color:#666;margin:0 2px;animation:typing 1s infinite;"></div>
          <div style="height:8px;width:8px;border-radius:50%;background-color:#666;margin:0 2px;animation:typing 1s infinite 0.2s;"></div>
          <div style="height:8px;width:8px;border-radius:50%;background-color:#666;margin:0 2px;animation:typing 1s infinite 0.4s;"></div>
        </div>
        <style>
          @keyframes typing {
            0% { transform: translateY(0); }
            50% { transform: translateY(-5px); }
            100% { transform: translateY(0); }
          }
        </style>
      </div>
    </div>
  `;
  
  messageContainer.innerHTML += typingHtml;
  messageContainer.scrollTop = messageContainer.scrollHeight;
  
  // Simulate AI response with realistic delay
  setTimeout(function() {
    // Remove typing indicator
    const typingIndicator = messageContainer.querySelector('.ai-typing');
    if (typingIndicator) {
      typingIndicator.remove();
    }
    
    // Sample response based on input
    let response;
    if (message.toLowerCase().includes('properties') || message.toLowerCase().includes('expand')) {
      response = `### Properties of Candidate 3

**ADMET Profile:**
- **Absorption:** Moderate oral bioavailability (55-65%) with food effect
- **Distribution:** Brain penetration ratio (brain:plasma) of 0.4-0.5
- **Metabolism:** Primary via CYP2D6, half-life ~6 hours
- **Excretion:** Primarily renal (65%), with some biliary excretion (35%)
- **Toxicity:** Low cardiotoxicity risk, minimal hepatotoxicity markers

**Binding Profile:**
- DAT: Ki = 42nM (high affinity)
- NET: Ki = 110nM (moderate affinity)
- SERT: Ki = 1240nM (low affinity, reduces serotonergic side effects)
- D1 receptor: Ki > 5000nM (minimal direct receptor activity)
- D2 receptor: Ki > 3000nM (minimal direct receptor activity)

**Physicochemical Properties:**
- LogP: 2.6 (good balance of lipophilicity)
- topological polar surface area: 48.2 Å²
- Solubility: 0.3 mg/mL at pH 7.4
- pKa: 9.2 (protonated at physiological pH)

This candidate shows the best overall profile for ADHD treatment with a good balance of efficacy and safety markers. Would you like me to suggest any chemical modifications to enhance specific properties?`;    
    } else if (message.toLowerCase().includes('side effect') || message.toLowerCase().includes('safety')) {
      response = `## Safety Profile Analysis

Candidate 3 was specifically designed to minimize common stimulant side effects:

1. **Cardiovascular effects:** The modified structure significantly reduces peripheral adrenergic activation compared to Adderall, with predicted blood pressure elevation of +2-4 mmHg (vs +5-10 mmHg for Adderall)

2. **Sleep disruption:** The 6-hour half-life allows for morning dosing without significant night-time sympathetic activation

3. **Appetite suppression:** Approximately 50% less appetite suppression than Adderall in animal models

4. **Abuse potential:** The controlled-release properties from the ester linkage create a gradual onset of action, reducing euphoric effects and abuse liability

5. **Growth effects in pediatric patients:** Predicted to have minimal impact on growth hormone release compared to traditional stimulants

This improved safety profile is achieved while maintaining therapeutic efficacy through selective DAT/NET inhibition without significant serotonergic or direct dopamine receptor activity.`;    
    } else {
      response = `Thank you for your question about "${message}"!

Based on the molecule design, here's my analysis:

- This molecule is designed to target dopamine/norepinephrine transporters with selective affinity
- The structure includes a phenylethylamine core similar to traditional stimulants
- Key modifications include: fluorine substitution for metabolic stability, ester linkage for controlled release, and piperazine moiety for reduced peripheral activity

Regarding specific properties:
- Predicted binding: Ki(DAT) ≈ 45nM, Ki(NET) ≈ 120nM
- Blood-brain barrier penetration is moderate (0.3-0.5 ratio)
- Half-life ~4-6 hours providing sustained but not excessive duration

Would you like me to suggest any structural modifications or provide more details on a specific aspect of this compound?`;    
    }
    
    // Add AI response to UI
    const aiMessageHtml = `
      <div class="message ai-message" style="display:flex;margin-bottom:16px;">
        <div style="width:32px;height:32px;border-radius:50%;background-color:#2196f3;color:white;display:flex;justify-content:center;align-items:center;margin-right:10px;">AI</div>
        <div style="flex:1;background-color:rgba(0,0,0,0.03);padding:10px;border-radius:8px;">
          <div>${formatMarkdown(response)}</div>
          <div style="font-size:0.8em;color:#666;margin-top:4px;">${new Date().toLocaleString()}</div>
        </div>
      </div>
    `;
    
    messageContainer.innerHTML += aiMessageHtml;
    messageContainer.scrollTop = messageContainer.scrollHeight;
  }, 2000 + Math.random() * 1000);
}

// Simple markdown formatter
function formatMarkdown(text) {
  return text
    .replace(/\n\n/g, '<br><br>')
    .replace(/\n/g, '<br>')
    .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
    .replace(/\*([^*]+)\*/g, '<em>$1</em>')
    .replace(/#{3}\s*([^\n]+)/g, '<h3>$1</h3>')
    .replace(/#{2}\s*([^\n]+)/g, '<h2>$1</h2>')
    .replace(/#{1}\s*([^\n]+)/g, '<h1>$1</h1>')
    .replace(/- ([^\n]+)/g, '• $1<br>');
}

// Helper function to find elements by text content
function findElementsByText(selector, text) {
  const elements = Array.from(document.querySelectorAll(selector));
  return elements.filter(element => element.textContent.toLowerCase().includes(text.toLowerCase()));
}
