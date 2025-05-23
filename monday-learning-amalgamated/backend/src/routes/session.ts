import { Router } from 'express'
import { logger } from '../utils/logger.js'

const router = Router()

// Session data store (in production, use Redis/database)
const sessions = new Map<string, any>()

// Start a new learning session
router.post('/start', async (req, res) => {
  try {
    const { userAgent, isVRSupported, timestamp } = req.body
    
    const sessionId = `session_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`
    
    const sessionData = {
      id: sessionId,
      startTime: timestamp || Date.now(),
      isActive: true,
      userAgent: userAgent || req.get('User-Agent'),
      isVRSupported: isVRSupported || false,
      totalQueries: 0,
      topicsExplored: [],
      learningMomentum: 0,
      hasCompletedIntro: false,
      lastInteractionTime: Date.now()
    }
    
    sessions.set(sessionId, sessionData)
    
    logger.info('Session started', { sessionId, userAgent: sessionData.userAgent })
    
    res.json({ 
      success: true, 
      message: 'Session started successfully',
      sessionId,
      data: sessionData
    })
  } catch (error) {
    logger.error('Failed to start session:', error)
    res.status(500).json({ success: false, error: 'Failed to start session' })
  }
})

// End a learning session
router.post('/end', async (req, res) => {
  try {
    const { sessionId } = req.body
    
    if (!sessionId || !sessions.has(sessionId)) {
      return res.status(404).json({ success: false, error: 'Session not found' })
    }
    
    const sessionData = sessions.get(sessionId)
    sessionData.isActive = false
    sessionData.endTime = Date.now()
    sessionData.duration = sessionData.endTime - sessionData.startTime
    
    sessions.set(sessionId, sessionData)
    
    logger.info('Session ended', { 
      sessionId, 
      duration: sessionData.duration,
      totalQueries: sessionData.totalQueries 
    })
    
    res.json({ 
      success: true, 
      message: 'Session ended successfully',
      data: sessionData
    })
  } catch (error) {
    logger.error('Failed to end session:', error)
    res.status(500).json({ success: false, error: 'Failed to end session' })
  }
})

// Get session status
router.get('/status/:sessionId', async (req, res) => {
  try {
    const { sessionId } = req.params
    
    if (!sessionId || !sessions.has(sessionId)) {
      return res.status(404).json({ success: false, error: 'Session not found' })
    }
    
    const sessionData = sessions.get(sessionId)
    
    res.json({ 
      success: true, 
      data: sessionData
    })
  } catch (error) {
    logger.error('Failed to get session status:', error)
    res.status(500).json({ success: false, error: 'Failed to get session status' })
  }
})

// Update session data
router.put('/update/:sessionId', async (req, res) => {
  try {
    const { sessionId } = req.params
    const updates = req.body
    
    if (!sessionId || !sessions.has(sessionId)) {
      return res.status(404).json({ success: false, error: 'Session not found' })
    }
    
    const sessionData = sessions.get(sessionId)
    Object.assign(sessionData, updates, { lastInteractionTime: Date.now() })
    sessions.set(sessionId, sessionData)
    
    res.json({ 
      success: true, 
      message: 'Session updated successfully',
      data: sessionData
    })
  } catch (error) {
    logger.error('Failed to update session:', error)
    res.status(500).json({ success: false, error: 'Failed to update session' })
  }
})

export default router 