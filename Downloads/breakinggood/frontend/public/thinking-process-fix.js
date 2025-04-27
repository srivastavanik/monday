// Original thinking process visualization with minimal fixes to prevent UI issues

document.addEventListener('DOMContentLoaded', function() {
  console.log('Applying minor UI fixes...');
  
  // Add critical CSS fixes
  const style = document.createElement('style');
  style.textContent = `
    /* Fix container height */
    .claude-is-designing {
      max-height: 100px;
      overflow: hidden;
    }
    
    /* Ensure molecules are always visible */
    .molecule-card {
      display: block !important;
      margin-bottom: 20px;
    }
    
    /* Prevent extremely long text from breaking layout */
    pre, code {
      max-height: 300px;
      overflow: auto;
    }
    
    /* Ensure chat containers remain the right size */
    .message-container {
      max-height: none !important;
    }
  `;
  
  document.head.appendChild(style);
  
  // Apply fixes at short intervals to ensure everything renders properly
  setInterval(function() {
    try {
      // Find all molecule cards and ensure they're visible
      const moleculeCards = document.querySelectorAll('.molecule-card');
      moleculeCards.forEach(card => {
        card.style.display = 'block';
        card.style.visibility = 'visible';
      });
      
      // Find chat container and make sure it's visible
      const chatContainer = document.querySelector('.message-container');
      if (chatContainer) {
        chatContainer.style.display = 'block';
        chatContainer.style.visibility = 'visible';
        chatContainer.style.maxHeight = 'none';
        chatContainer.style.overflow = 'visible';
      }
      
      // Find result container and ensure it's visible
      const resultContainer = document.querySelector('.generated-molecules');
      if (resultContainer) {
        resultContainer.style.display = 'block';
        resultContainer.style.visibility = 'visible';
        resultContainer.style.maxHeight = 'none';
      }
    } catch (err) {
      console.error('Error in display fix:', err);
    }
  }, 200);
});
