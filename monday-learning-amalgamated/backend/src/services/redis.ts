import { getRedisClient } from '../database/index.js'
import { logger } from '../utils/logger.js'

export async function initializeRedis() {
  try {
    const client = getRedisClient()
    if (!client) {
      throw new Error('Redis client not initialized')
    }
    logger.info('Redis service initialized')
  } catch (error) {
    logger.error('Redis service initialization failed:', error)
    throw error
  }
}

// Session management
export class SessionManager {
  private static getSessionKey(sessionId: string): string {
    return `session:${sessionId}`
  }

  static async createSession(sessionId: string, data: any): Promise<void> {
    try {
      const client = getRedisClient()
      if (!client) throw new Error('Redis client not available')
      
      const sessionData = {
        ...data,
        createdAt: Date.now(),
        lastAccessed: Date.now()
      }
      
      await client.setEx(
        this.getSessionKey(sessionId),
        3600, // 1 hour TTL
        JSON.stringify(sessionData)
      )
      
      logger.info('Session created in Redis', { sessionId })
    } catch (error) {
      logger.error('Failed to create session in Redis:', error)
      throw error
    }
  }

  static async getSession(sessionId: string): Promise<any> {
    try {
      const client = getRedisClient()
      if (!client) return null
      
      const data = await client.get(this.getSessionKey(sessionId))
      if (!data) return null
      
      const sessionData = JSON.parse(data)
      
      // Update last accessed time
      sessionData.lastAccessed = Date.now()
      await client.setEx(
        this.getSessionKey(sessionId),
        3600,
        JSON.stringify(sessionData)
      )
      
      return sessionData
    } catch (error) {
      logger.error('Failed to get session from Redis:', error)
      return null
    }
  }

  static async updateSession(sessionId: string, updates: any): Promise<void> {
    try {
      const client = getRedisClient()
      if (!client) throw new Error('Redis client not available')
      
      const existingData = await this.getSession(sessionId)
      if (!existingData) {
        throw new Error('Session not found')
      }
      
      const updatedData = {
        ...existingData,
        ...updates,
        lastAccessed: Date.now()
      }
      
      await client.setEx(
        this.getSessionKey(sessionId),
        3600,
        JSON.stringify(updatedData)
      )
      
      logger.info('Session updated in Redis', { sessionId })
    } catch (error) {
      logger.error('Failed to update session in Redis:', error)
      throw error
    }
  }

  static async deleteSession(sessionId: string): Promise<void> {
    try {
      const client = getRedisClient()
      if (!client) return
      
      await client.del(this.getSessionKey(sessionId))
      logger.info('Session deleted from Redis', { sessionId })
    } catch (error) {
      logger.error('Failed to delete session from Redis:', error)
      throw error
    }
  }
}

// Cache management
export class CacheManager {
  static async set(key: string, value: any, ttl: number = 3600): Promise<void> {
    try {
      const client = getRedisClient()
      if (!client) return
      
      await client.setEx(key, ttl, JSON.stringify(value))
    } catch (error) {
      logger.error('Failed to set cache value:', error)
    }
  }

  static async get(key: string): Promise<any> {
    try {
      const client = getRedisClient()
      if (!client) return null
      
      const data = await client.get(key)
      return data ? JSON.parse(data) : null
    } catch (error) {
      logger.error('Failed to get cache value:', error)
      return null
    }
  }

  static async delete(key: string): Promise<void> {
    try {
      const client = getRedisClient()
      if (!client) return
      
      await client.del(key)
    } catch (error) {
      logger.error('Failed to delete cache value:', error)
    }
  }

  static async exists(key: string): Promise<boolean> {
    try {
      const client = getRedisClient()
      if (!client) return false
      
      const result = await client.exists(key)
      return result === 1
    } catch (error) {
      logger.error('Failed to check cache existence:', error)
      return false
    }
  }
} 