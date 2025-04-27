// React DOM manipulation safety fix
// This script patches DOM methods to prevent common React errors

console.log('Applying React DOM safety patches...');

(function() {
  // Store original DOM methods
  const originalRemoveChild = Node.prototype.removeChild;
  const originalInsertBefore = Node.prototype.insertBefore;
  const originalAppendChild = Node.prototype.appendChild;
  const originalReplaceChild = Node.prototype.replaceChild;
  
  // Patch removeChild to handle errors gracefully
  Node.prototype.removeChild = function(child) {
    try {
      // Check if child is actually a child of this node
      if (!this.contains(child)) {
        console.warn('React DOM fix: Attempted to remove a child that is not part of the parent');
        return child; // Return the child without throwing an error
      }
      return originalRemoveChild.call(this, child);
    } catch (e) {
      console.warn('React DOM fix: Handled error in removeChild', e);
      return child;
    }
  };
  
  // Patch insertBefore to handle errors gracefully
  Node.prototype.insertBefore = function(newNode, referenceNode) {
    try {
      // If reference node doesn't exist or isn't a child, append instead
      if (referenceNode && !this.contains(referenceNode)) {
        console.warn('React DOM fix: Reference node not found in insertBefore, appending instead');
        return originalAppendChild.call(this, newNode);
      }
      return originalInsertBefore.call(this, newNode, referenceNode);
    } catch (e) {
      console.warn('React DOM fix: Handled error in insertBefore', e);
      try {
        return originalAppendChild.call(this, newNode);
      } catch (appendError) {
        console.warn('React DOM fix: Also failed to append as fallback', appendError);
        return newNode;
      }
    }
  };
  
  // Patch replaceChild to handle errors gracefully
  Node.prototype.replaceChild = function(newChild, oldChild) {
    try {
      // Check if oldChild is actually a child of this node
      if (!this.contains(oldChild)) {
        console.warn('React DOM fix: Attempted to replace a child that is not part of the parent');
        try {
          return originalAppendChild.call(this, newChild);
        } catch (appendError) {
          console.warn('React DOM fix: Failed to append as fallback', appendError);
          return newChild;
        }
      }
      return originalReplaceChild.call(this, newChild, oldChild);
    } catch (e) {
      console.warn('React DOM fix: Handled error in replaceChild', e);
      try {
        return originalAppendChild.call(this, newChild);
      } catch (appendError) {
        console.warn('React DOM fix: Also failed to append as fallback', appendError);
        return newChild;
      }
    }
  };
  
  console.log('React DOM safety patches applied successfully');
})();
