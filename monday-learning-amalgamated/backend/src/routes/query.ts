import { Router } from 'express'
import { perplexityService } from '../services/perplexity.js'
import { logger } from '../utils/logger.js'

const router = Router()

// Basic query endpoint
router.post('/basic', async (req, res) => {
  try {
    const { query, context, sessionId } = req.body
    
    if (!query) {
      return res.status(400).json({ success: false, error: 'Query is required' })
    }
    
    const response = await perplexityService.processQuery({
      query,
      mode: 'basic',
      context,
      sessionId
    })
    
    logger.info('Basic query processed', { sessionId, query: query.substring(0, 50) })
    
    res.json({ 
      success: true, 
      data: response,
      message: 'Basic query processed successfully'
    })
  } catch (error) {
    logger.error('Basic query failed:', error)
    res.status(500).json({ success: false, error: 'Failed to process basic query' })
  }
})

// Reasoning query endpoint
router.post('/reasoning', async (req, res) => {
  try {
    const { query, context, sessionId } = req.body
    
    if (!query) {
      return res.status(400).json({ success: false, error: 'Query is required' })
    }
    
    const response = await perplexityService.processQuery({
      query,
      mode: 'reasoning',
      context,
      sessionId
    })
    
    logger.info('Reasoning query processed', { sessionId, query: query.substring(0, 50) })
    
    res.json({ 
      success: true, 
      data: response,
      message: 'Reasoning query processed successfully'
    })
  } catch (error) {
    logger.error('Reasoning query failed:', error)
    res.status(500).json({ success: false, error: 'Failed to process reasoning query' })
  }
})

// Deep research query endpoint
router.post('/research', async (req, res) => {
  try {
    const { query, sessionId } = req.body
    
    if (!query) {
      return res.status(400).json({ success: false, error: 'Query is required' })
    }
    
    const response = await perplexityService.processQuery({
      query,
      mode: 'research',
      sessionId
    })
    
    logger.info('Research query processed', { sessionId, query: query.substring(0, 50) })
    
    res.json({ 
      success: true, 
      data: response,
      message: 'Research query processed successfully'
    })
  } catch (error) {
    logger.error('Research query failed:', error)
    res.status(500).json({ success: false, error: 'Failed to process research query' })
  }
})

// Get query history for a session
router.get('/history/:sessionId', async (req, res) => {
  try {
    const { sessionId } = req.params
    const { limit = 10 } = req.query
    
    // In production, fetch from database
    const history = []
    
    res.json({ 
      success: true, 
      data: {
        sessionId,
        queries: history,
        total: history.length
      }
    })
  } catch (error) {
    logger.error('Failed to get query history:', error)
    res.status(500).json({ success: false, error: 'Failed to get query history' })
  }
})

export default router 