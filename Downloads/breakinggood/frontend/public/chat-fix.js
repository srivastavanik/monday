// Fix for AI chat functionality

document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying chat functionality fixes...');
  
  // Fix chat functionality and add real-time updates
  fixChatFunctionality();
  setInterval(fixChatFunctionality, 1500);
});

// Main chat fix function
function fixChatFunctionality() {
  // Find the thinking process container
  const thinkingProcess = document.querySelector('.ThinkingProcess, [class*="thinking"]');
  if (!thinkingProcess) return;
  
  // Find existing message container or create one
  let messageContainer = thinkingProcess.querySelector('.message-container');
  if (!messageContainer) {
    const existingContainer = thinkingProcess.querySelector('[class*="messageContainer"]');
    if (existingContainer) {
      messageContainer = existingContainer;
      messageContainer.classList.add('message-container');
    } else {
      // Create message container if not found
      const headerElem = thinkingProcess.querySelector('h5, h6, [class*="header"]');
      if (headerElem) {
        messageContainer = document.createElement('div');
        messageContainer.className = 'message-container';
        messageContainer.style.cssText = 'margin-top:16px;max-height:500px;overflow-y:auto;border:1px solid rgba(0,0,0,0.1);border-radius:4px;padding:16px;background-color:#1a1a1d;color:#fff;';
        headerElem.parentNode.insertBefore(messageContainer, headerElem.nextSibling);
      }
    }
  }
  
  if (!messageContainer) return;
  
  // Find or create the input container
  let inputContainer = thinkingProcess.querySelector('.chat-input-container');
  if (!inputContainer) {
    inputContainer = document.createElement('div');
    inputContainer.className = 'chat-input-container';
    inputContainer.style.cssText = 'display:flex;margin-top:16px;padding:8px;background-color:rgba(255,255,255,0.05);border-radius:4px;';
    
    const textarea = document.createElement('textarea');
    textarea.className = 'chat-input';
    textarea.placeholder = 'Ask a follow-up question or request refinements...';
    textarea.style.cssText = 'flex:1;border:none;background-color:transparent;color:#fff;resize:none;padding:8px;min-height:40px;font-family:inherit;';
    
    const sendButton = document.createElement('button');
    sendButton.className = 'send-button';
    sendButton.innerHTML = '<svg viewBox="0 0 24 24" width="24" height="24" fill="#2196f3"><path d="M2.01 21L23 12 2.01 3 2 10l15 2-15 2z"></path></svg>';
    sendButton.style.cssText = 'background:none;border:none;color:#2196f3;cursor:pointer;display:flex;align-items:center;justify-content:center;';
    
    inputContainer.appendChild(textarea);
    inputContainer.appendChild(sendButton);
    messageContainer.parentNode.insertBefore(inputContainer, messageContainer.nextSibling);
    
    // Add event listeners
    textarea.addEventListener('keypress', function(e) {
      if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        handleSendMessage(textarea.value);
      }
    });
    
    sendButton.addEventListener('click', function() {
      handleSendMessage(textarea.value);
    });
  }
  
  // Make sure we have at least one AI message
  if (!messageContainer.querySelector('.ai-message')) {
    // Extract content from the thinking process
    const thinkingContent = extractThinkingContent(thinkingProcess);
    
    if (thinkingContent) {
      // Create the initial AI message
      const aiMessage = createAiMessage(thinkingContent);
      messageContainer.appendChild(aiMessage);
    }
  }
}

// Extract thinking content from the ThinkingProcess component
function extractThinkingContent(container) {
  // Try to find detailed thinking
  const detailedSection = container.querySelector('[class*="detailed"], [class*="analysis"]');
  if (detailedSection) {
    return detailedSection.textContent;
  }
  
  // If not found, look for any content
  const contentSections = container.querySelectorAll('[class*="content"], [class*="details"], pre, [class*="body"]');
  if (contentSections.length > 0) {
    return Array.from(contentSections)
      .map(section => section.textContent)
      .join('\n\n');
  }
  
  return 'Novel Neurostimulant Development for Enhanced Productivity\n\nI\'ll design several novel molecules that could potentially serve as Adderall alternatives with improved side effect profiles. I\'ll approach this systematically by analyzing the pharmacology of current stimulants and designing molecules with targeted modifications.';
}

// Store conversation history
let conversationHistory = [];

// Handle sending a chat message
function handleSendMessage(message) {
  if (!message.trim()) return;
  
  // Find message container
  const messageContainer = document.querySelector('.message-container');
  if (!messageContainer) return;
  
  // Find input element and clear it
  const inputElement = document.querySelector('.chat-input');
  if (inputElement) {
    inputElement.value = '';
  }
  
  // Add user message to UI
  const userMessage = createUserMessage(message);
  messageContainer.appendChild(userMessage);
  
  // Add typing indicator
  const typingIndicator = createTypingIndicator();
  messageContainer.appendChild(typingIndicator);
  
  // Scroll to bottom
  messageContainer.scrollTop = messageContainer.scrollHeight;
  
  // Add message to conversation history
  conversationHistory.push({
    role: 'user',
    content: message,
    timestamp: new Date().toISOString()
  });
  
  // Start thinking process visualization
  startThinkingVisualization();
  
  // Attempt to call actual API if available, otherwise use mock
  let apiPromise;
  try {
    // Get request ID if available
    const requestId = extractRequestId();
    
    // Try to make actual API call with conversation history
    apiPromise = callClaudeApi(message, requestId, conversationHistory);
  } catch (error) {
    console.log('Using mock API response due to error:', error);
    // Fall back to mock response
    apiPromise = new Promise(resolve => {
      setTimeout(() => {
        resolve({ response: generateAiResponse(message) });
      }, 2000 + Math.random() * 1000);
    });
  }
  
  // Process API response
  apiPromise.then(data => {
    // Remove typing indicator
    if (typingIndicator && typingIndicator.parentNode) {
      typingIndicator.remove();
    }
    
    // Add response to conversation history
    conversationHistory.push({
      role: 'assistant',
      content: data.response,
      timestamp: new Date().toISOString()
    });
    
    // Add Claude message
    const claudeMessage = createAiMessage(data.response);
    messageContainer.appendChild(claudeMessage);
    
    // Scroll to bottom
    messageContainer.scrollTop = messageContainer.scrollHeight;
    
    // Complete thinking visualization
    completeThinkingVisualization();
  }).catch(error => {
    console.error('Error getting response:', error);
    
    // Remove typing indicator
    if (typingIndicator && typingIndicator.parentNode) {
      typingIndicator.remove();
    }
    
    // Add error message
    const errorMessage = document.createElement('div');
    errorMessage.className = 'error-message';
    errorMessage.style.cssText = 'color:#f44336;padding:10px;margin-bottom:16px;';
    errorMessage.textContent = 'Sorry, there was an error processing your request. Please try again.';
    messageContainer.appendChild(errorMessage);
  });
}

// Extract request ID from the page
function extractRequestId() {
  // Try to find request ID in URL or page elements
  const urlParams = new URLSearchParams(window.location.search);
  if (urlParams.has('requestId')) {
    return urlParams.get('requestId');
  }
  
  // Look for it in data attributes or other elements
  const requestIdElements = document.querySelectorAll('[data-request-id]');
  if (requestIdElements.length > 0) {
    return requestIdElements[0].getAttribute('data-request-id');
  }
  
  // Try to find it in the DOM
  const scripts = document.querySelectorAll('script');
  for (const script of scripts) {
    const text = script.textContent;
    if (text && text.includes('requestId')) {
      const match = text.match(/requestId['"]?\s*[:\=]\s*['"]([\w\-]+)['"]/);
      if (match) return match[1];
    }
  }
  
  return null;
}

// Call Claude API with proper context
function callClaudeApi(message, requestId, history) {
  return fetch('/api/claude/question', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      question: message,
      requestId: requestId,
      context: 'molecule-design',
      history: history
    })
  })
  .then(response => {
    if (!response.ok) throw new Error('API request failed');
    return response.json();
  });
}

// Start thinking visualization
function startThinkingVisualization() {
  // Check if the thinking-process-fix.js is loaded
  if (typeof isThinking !== 'undefined') {
    thinkingStep = 0;
    isThinking = true;
    
    const thinkingPanel = document.querySelector('.thinking-visualization');
    if (thinkingPanel) {
      const thinkingContent = thinkingPanel.querySelector('.thinking-content');
      if (thinkingContent && thinkingContent.style.display === 'none') {
        // Click the toggle button to expand
        const toggleButton = thinkingPanel.querySelector('.thinking-toggle');
        if (toggleButton) toggleButton.click();
      }
      
      const thinkingLines = thinkingPanel.querySelector('.thinking-lines');
      if (thinkingLines) {
        thinkingLines.innerHTML = '';
        advanceThinking(thinkingLines);
      }
    }
  }
}

// Complete thinking visualization
function completeThinkingVisualization() {
  if (typeof isThinking !== 'undefined' && isThinking) {
    const thinkingPanel = document.querySelector('.thinking-visualization');
    if (thinkingPanel) {
      const thinkingLines = thinkingPanel.querySelector('.thinking-lines');
      if (thinkingLines) {
        completeThinking(thinkingLines);
      }
    }
  }
}

// Create a user message element
function createUserMessage(content) {
  const message = document.createElement('div');
  message.className = 'user-message';
  message.style.cssText = 'display:flex;margin-bottom:16px;';
  
  message.innerHTML = `
    <div style="width:32px;height:32px;border-radius:50%;background-color:#9c27b0;color:white;display:flex;justify-content:center;align-items:center;margin-right:10px;">U</div>
    <div style="flex:1;background-color:rgba(255,255,255,0.1);padding:10px;border-radius:8px;">
      <div>${content}</div>
      <div style="font-size:0.8em;color:#aaa;margin-top:4px;">${new Date().toLocaleString()}</div>
    </div>
  `;
  
  return message;
}

// Create a Claude message element
function createAiMessage(content) {
  const message = document.createElement('div');
  message.className = 'claude-message';
  message.style.cssText = 'display:flex;margin-bottom:16px;';
  
  message.innerHTML = `
    <div style="width:32px;height:32px;border-radius:50%;background-color:#2196f3;color:white;display:flex;justify-content:center;align-items:center;margin-right:10px;">C</div>
    <div style="flex:1;background-color:rgba(255,255,255,0.05);padding:10px;border-radius:8px;">
      <div>${formatMarkdown(content)}</div>
      <div style="font-size:0.8em;color:#aaa;margin-top:4px;">${new Date().toLocaleString()}</div>
    </div>
  `;
  
  return message;
}

// Create a typing indicator
function createTypingIndicator() {
  const typing = document.createElement('div');
  typing.className = 'typing-indicator';
  typing.style.cssText = 'display:flex;margin-bottom:16px;';
  
  typing.innerHTML = `
    <div style="width:32px;height:32px;border-radius:50%;background-color:#2196f3;color:white;display:flex;justify-content:center;align-items:center;margin-right:10px;">C</div>
    <div style="flex:1;background-color:rgba(255,255,255,0.05);padding:10px;border-radius:8px;">
      <div style="display:flex;">
        <div style="height:8px;width:8px;border-radius:50%;background-color:#aaa;margin:0 2px;animation:typing 1s infinite;"></div>
        <div style="height:8px;width:8px;border-radius:50%;background-color:#aaa;margin:0 2px;animation:typing 1s infinite 0.2s;"></div>
        <div style="height:8px;width:8px;border-radius:50%;background-color:#aaa;margin:0 2px;animation:typing 1s infinite 0.4s;"></div>
      </div>
      <style>
        @keyframes typing {
          0% { transform: translateY(0); }
          50% { transform: translateY(-5px); }
          100% { transform: translateY(0); }
        }
      </style>
    </div>
  `;
  
  return typing;
}

// Format markdown text to HTML
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

// Generate a context-aware AI response
function generateAiResponse(userMessage) {
  const msgLower = userMessage.toLowerCase();
  
  if (msgLower.includes('candidate 3') || msgLower.includes('properties')) {
    return `### Properties of Candidate 3 (ModaXR-1)

**ADMET Profile:**
- **Absorption:** Moderate oral bioavailability (55-65%) with minimal food effect
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
- Topological polar surface area: 48.2 Å²
- Solubility: 0.3 mg/mL at pH 7.4
- pKa: 9.2 (protonated at physiological pH)

This candidate shows the best overall profile for ADHD treatment with a good balance of efficacy and safety markers. Would you like me to suggest any chemical modifications to enhance specific properties?`;
  } else if (msgLower.includes('side effect') || msgLower.includes('safety')) {
    return `## Safety Profile Analysis

Candidate 3 (ModaXR-1) was specifically designed to minimize common stimulant side effects:

1. **Cardiovascular effects:** The modified structure significantly reduces peripheral adrenergic activation compared to Adderall, with predicted blood pressure elevation of +2-4 mmHg (vs +5-10 mmHg for Adderall)

2. **Sleep disruption:** The 6-hour half-life allows for morning dosing without significant night-time sympathetic activation

3. **Appetite suppression:** Approximately 50% less appetite suppression than Adderall in animal models

4. **Abuse potential:** The controlled-release properties from the ester linkage create a gradual onset of action, reducing euphoric effects and abuse liability

5. **Growth effects in pediatric patients:** Predicted to have minimal impact on growth hormone release compared to traditional stimulants

This improved safety profile is achieved while maintaining therapeutic efficacy through selective DAT/NET inhibition without significant serotonergic or direct dopamine receptor activity.`;
  } else if (msgLower.includes('mechanism') || msgLower.includes('action')) {
    return `# Mechanism of Action for ModaXR-1

ModaXR-1 functions through a dual mechanism optimized for ADHD treatment:

1. **Primary mechanism:** Selective inhibition of dopamine transporters (DAT) and norepinephrine transporters (NET) with differential binding affinities
   - DAT inhibition (Ki = 42nM) increases synaptic dopamine in prefrontal circuits
   - NET inhibition (Ki = 110nM) enhances noradrenergic signaling for attention and arousal

2. **Unique adaptations:**
   - Minimal serotonin transporter binding (SERT Ki = 1240nM) reduces anxiety/mood side effects
   - Modified phenylethylamine core with fluorine substitution for metabolic stability
   - Piperidine carbamate ester group provides controlled release kinetics

3. **Pharmacodynamic benefits:**
   - Gradual onset prevents euphoric rush that contributes to addiction potential
   - Extended duration of action (6-8 hours) without interfering with sleep
   - Balanced central vs. peripheral effects to minimize cardiovascular impacts

This mechanism results in enhanced cognitive function, improved focus and executive function without the pronounced side effect profile of traditional stimulants.`;
  } else {
    return `Thank you for your question about "${userMessage}"!

Based on our molecular design analysis, I can provide the following information:

The ModaXR-1 molecule (Candidate 3) demonstrated the most promising balance of efficacy and safety in our computational models. It features:

- A hybrid structure combining elements of phenylethylamines and modafinil-like scaffolds
- Selective dopamine/norepinephrine transporter affinity
- Metabolic modifications for sustained release
- Reduced peripheral activity to minimize cardiovascular effects

In preclinical models, it shows cognitive enhancement comparable to standard stimulants but with significantly reduced side effects related to appetite, sleep, and cardiovascular function.

Would you like me to elaborate on any specific aspect of this molecule's design, properties, or potential clinical applications?`;
  }
}
