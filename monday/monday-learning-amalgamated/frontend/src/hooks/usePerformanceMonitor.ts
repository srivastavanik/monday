import { useState, useEffect, useRef } from 'react'

interface PerformanceMetrics {
  fps: number
  frameTime: number
  memoryUsage: number
}

export const usePerformanceMonitor = (): PerformanceMetrics => {
  const [fps, setFps] = useState(60)
  const [frameTime, setFrameTime] = useState(16.67)
  const [memoryUsage, setMemoryUsage] = useState(0)
  
  const frameCountRef = useRef(0)
  const lastTimeRef = useRef(performance.now())

  useEffect(() => {
    const measurePerformance = () => {
      const now = performance.now()
      const delta = now - lastTimeRef.current
      
      frameCountRef.current++
      
      if (delta >= 1000) {
        const currentFps = Math.round((frameCountRef.current * 1000) / delta)
        const currentFrameTime = delta / frameCountRef.current
        
        setFps(currentFps)
        setFrameTime(currentFrameTime)
        
        frameCountRef.current = 0
        lastTimeRef.current = now
      }
      
      if ('memory' in performance) {
        const memory = (performance as any).memory
        setMemoryUsage(memory.usedJSHeapSize / 1024 / 1024)
      }
      
      requestAnimationFrame(measurePerformance)
    }

    requestAnimationFrame(measurePerformance)
  }, [])

  return { fps, frameTime, memoryUsage }
} 