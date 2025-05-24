import { Router } from 'express'
import { logger } from '../utils/logger.js'

// Define interfaces for analytics data structures
interface Session {
  id: string;
  startTime: number;
  endTime?: number;
  // Add other session-specific fields if needed
}

interface Query {
  id: string;
  sessionId: string;
  timestamp: number;
  type: 'basic' | 'reasoning' | 'research'; // Example query types
  data: any; // This can be more specific based on actual query data structure
  responseTime?: number;
}

interface PerformanceEntry {
  id: string;
  timestamp: number;
  fps: number;
  frameTime: number;
  memoryUsage: number;
  // Add other performance metrics as needed
}

interface Interaction {
  id: string;
  sessionId: string;
  timestamp: number;
  type: string; // e.g., 'voice_command', 'panel_manipulation'
  data: any; // Specific data for the interaction type
}

interface AnalyticsStore {
  sessions: Session[];
  queries: Query[];
  performance: PerformanceEntry[];
  interactions: Interaction[];
}

const router = Router()

// In-memory analytics store (in production, use database)
const analytics: AnalyticsStore = {
  sessions: [],
  queries: [],
  performance: [],
  interactions: []
}

// Get analytics metrics
router.get('/metrics', async (req, res) => {
  try {
    const metrics = {
      totalSessions: analytics.sessions.length,
      totalQueries: analytics.queries.length,
      averagePerformance: {
        fps: 72,
        frameTime: 13.9,
        memoryUsage: 256
      }
    }
    
    res.json({ 
      success: true, 
      data: metrics
    })
  } catch (error) {
    logger.error('Failed to get analytics metrics:', error)
    res.status(500).json({ success: false, error: 'Failed to get analytics metrics' })
  }
})

// Track an event
router.post('/track', async (req, res) => {
  try {
    const { type, data, sessionId } = req.body
    
    const event = {
      id: `event_${Date.now()}`,
      type,
      data,
      sessionId,
      timestamp: Date.now()
    }
    
    analytics.interactions.push(event)
    
    res.json({ 
      success: true, 
      message: 'Event tracked successfully',
      eventId: event.id
    })
  } catch (error) {
    logger.error('Failed to track event:', error)
    res.status(500).json({ success: false, error: 'Failed to track event' })
  }
})

// Get learning insights for a session
router.get('/insights/:sessionId', async (req, res) => {
  try {
    const { sessionId } = req.params
    
    const sessionQueries = analytics.queries.filter(q => q.sessionId === sessionId)
    const sessionInteractions = analytics.interactions.filter(i => i.sessionId === sessionId)
    
    const insights = {
      totalQueries: sessionQueries.length,
      topicsExplored: [...new Set(sessionQueries.map(q => q.data.topic).filter(Boolean))],
      queryTypes: {
        basic: sessionQueries.filter(q => q.data.mode === 'basic').length,
        reasoning: sessionQueries.filter(q => q.data.mode === 'reasoning').length,
        research: sessionQueries.filter(q => q.data.mode === 'research').length
      },
      learningPath: sessionQueries.map(q => ({
        timestamp: q.timestamp,
        topic: q.data.topic,
        mode: q.data.mode,
        duration: q.data.responseTime
      })),
      engagementScore: calculateEngagementScore(sessionQueries, sessionInteractions)
    }
    
    res.json({ 
      success: true, 
      data: insights,
      sessionId
    })
  } catch (error) {
    logger.error('Failed to get learning insights:', error)
    res.status(500).json({ success: false, error: 'Failed to get learning insights' })
  }
})

// Helper functions
function getTopTopics(timeFilter: number) {
  const recentQueries = analytics.queries.filter(q => q.timestamp > timeFilter)
  const topicCounts = {}
  
  recentQueries.forEach(q => {
    const topic = q.data?.topic
    if (topic) {
      topicCounts[topic] = (topicCounts[topic] || 0) + 1
    }
  })
  
  return Object.entries(topicCounts)
    .sort(([,a], [,b]) => (b as number) - (a as number))
    .slice(0, 10)
    .map(([topic, count]) => ({ topic, count }))
}

function getLearningProgression(timeFilter: number) {
  const recentQueries = analytics.queries.filter(q => q.timestamp > timeFilter)
  const hourlyData = {}
  
  recentQueries.forEach(q => {
    const hour = Math.floor(q.timestamp / (60 * 60 * 1000)) * (60 * 60 * 1000)
    hourlyData[hour] = (hourlyData[hour] || 0) + 1
  })
  
  return Object.entries(hourlyData)
    .sort(([a], [b]) => Number(a) - Number(b))
    .map(([timestamp, count]) => ({ timestamp: Number(timestamp), queries: count }))
}

function calculateEngagementScore(queries: any[], interactions: any[]) {
  // Simple engagement score based on query frequency, variety, and interaction types
  const queryVariety = new Set(queries.map(q => q.data?.mode)).size
  const interactionVariety = new Set(interactions.map(i => i.type)).size
  const frequency = queries.length
  
  return Math.min(100, (queryVariety * 20) + (interactionVariety * 15) + (frequency * 2))
}

export default router 