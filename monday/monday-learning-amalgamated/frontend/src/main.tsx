import React from 'react'
import ReactDOM from 'react-dom/client'
import App from './App'
import './index.css'
import 'aframe'

// WebXR polyfill detection and setup
const initializeApp = async () => {
  // Check for WebXR support
  if (!navigator.xr) {
    console.warn('WebXR not supported, running in fallback mode')
  } else {
    try {
      const isSupported = await navigator.xr.isSessionSupported('immersive-vr')
      if (!isSupported) {
        console.warn('Immersive VR not supported, running in fallback mode')
      }
    } catch (error) {
      console.warn('WebXR session check failed:', error)
    }
  }

  // Remove loading screen
  const loadingElement = document.getElementById('loading')
  if (loadingElement) {
    loadingElement.style.opacity = '0'
    setTimeout(() => {
      loadingElement.remove()
    }, 500)
  }

  // Mount React app
  const root = ReactDOM.createRoot(document.getElementById('root')!)
  root.render(
    <React.StrictMode>
      <App />
    </React.StrictMode>
  )
}

// Initialize when DOM is ready
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', initializeApp)
} else {
  initializeApp()
} 