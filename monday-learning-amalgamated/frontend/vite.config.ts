import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { VitePWA } from 'vite-plugin-pwa'

export default defineConfig({
  plugins: [
    react(),
    VitePWA({
      registerType: 'autoUpdate',
      workbox: {
        globPatterns: ['**/*.{js,css,html,ico,png,svg}']
      },
      manifest: {
        name: 'Monday: Learning, Amalgamated',
        short_name: 'Monday',
        description: 'Voice-driven AI learning companion for VR',
        theme_color: '#20808D',
        background_color: '#091717',
        display: 'standalone',
        icons: [
          {
            src: 'icon-192x192.png',
            sizes: '192x192',
            type: 'image/png'
          },
          {
            src: 'icon-512x512.png',
            sizes: '512x512',
            type: 'image/png'
          }
        ]
      }
    })
  ],
  server: {
    host: '0.0.0.0',
    port: 3000,
    https: {
      // Self-signed certificate for local development
      // Quest requires HTTPS for WebXR and microphone access
      key: './certs/key.pem',
      cert: './certs/cert.pem'
    },
    proxy: {
      '/api': {
        target: 'http://localhost:3001',
        changeOrigin: true
      },
      '/socket.io': {
        target: 'http://localhost:3001',
        ws: true
      }
    }
  },
  build: {
    target: 'es2020',
    rollupOptions: {
      output: {
        manualChunks: {
          vendor: ['react', 'react-dom'],
          three: ['three', '@react-three/fiber', '@react-three/drei'],
          aframe: ['aframe', 'aframe-react'],
          audio: ['tone', 'wavesurfer.js']
        }
      }
    }
  },
  optimizeDeps: {
    include: ['aframe', 'three']
  },
  define: {
    global: 'globalThis'
  }
}) 