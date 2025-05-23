import { Client as PostgresClient, Pool } from 'pg'
import { createClient as createRedisClient } from 'redis'
import neo4j from 'neo4j-driver'
import { logger } from '../utils/logger.js'

let postgresPool: Pool | null = null
let redisClient: any = null
let neo4jDriver: any = null

// PostgreSQL initialization
async function initializePostgreSQL() {
  try {
    postgresPool = new Pool({
      connectionString: process.env.DATABASE_URL,
      ssl: process.env.NODE_ENV === 'production' ? { rejectUnauthorized: false } : false,
      max: 20,
      idleTimeoutMillis: 30000,
      connectionTimeoutMillis: 2000,
    })

    // Test connection
    const client = await postgresPool.connect()
    await client.query('SELECT NOW()')
    client.release()

    // Create tables if they don't exist
    await createTables()
    
    logger.info('PostgreSQL connected successfully')
  } catch (error) {
    logger.error('PostgreSQL connection failed:', error)
    throw error
  }
}

// Redis initialization
async function initializeRedis() {
  try {
    redisClient = createRedisClient({
      url: process.env.REDIS_URL
    })

    redisClient.on('error', (err: any) => {
      logger.error('Redis client error:', err)
    })

    redisClient.on('connect', () => {
      logger.info('Redis connected successfully')
    })

    await redisClient.connect()
  } catch (error) {
    logger.error('Redis connection failed:', error)
    throw error
  }
}

// Neo4j initialization
async function initializeNeo4j() {
  try {
    neo4jDriver = neo4j.driver(
      process.env.NEO4J_URI || 'bolt://localhost:7687',
      neo4j.auth.basic(
        process.env.NEO4J_USER || 'neo4j',
        process.env.NEO4J_PASSWORD || 'password'
      )
    )

    // Test connection
    const session = neo4jDriver.session()
    await session.run('RETURN 1 as test')
    await session.close()

    logger.info('Neo4j connected successfully')
  } catch (error) {
    logger.error('Neo4j connection failed:', error)
    throw error
  }
}

// Create PostgreSQL tables
async function createTables() {
  if (!postgresPool) return

  const queries = [
    // Users table
    `CREATE TABLE IF NOT EXISTS users (
      id SERIAL PRIMARY KEY,
      username VARCHAR(255) UNIQUE,
      email VARCHAR(255) UNIQUE,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      last_login TIMESTAMP
    )`,

    // Sessions table
    `CREATE TABLE IF NOT EXISTS sessions (
      id VARCHAR(255) PRIMARY KEY,
      user_id INTEGER REFERENCES users(id),
      start_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      end_time TIMESTAMP,
      is_active BOOLEAN DEFAULT TRUE,
      user_agent TEXT,
      is_vr_supported BOOLEAN DEFAULT FALSE,
      total_queries INTEGER DEFAULT 0,
      learning_momentum INTEGER DEFAULT 0,
      topics_explored TEXT[],
      session_data JSONB
    )`,

    // Queries table
    `CREATE TABLE IF NOT EXISTS queries (
      id SERIAL PRIMARY KEY,
      session_id VARCHAR(255) REFERENCES sessions(id),
      query_text TEXT NOT NULL,
      query_mode VARCHAR(50) NOT NULL,
      response_content TEXT,
      response_time_ms INTEGER,
      tokens_used INTEGER,
      citations JSONB,
      reasoning_steps JSONB,
      created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )`,

    // Analytics table
    `CREATE TABLE IF NOT EXISTS analytics (
      id SERIAL PRIMARY KEY,
      session_id VARCHAR(255) REFERENCES sessions(id),
      event_type VARCHAR(100) NOT NULL,
      event_data JSONB,
      timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
      user_agent TEXT,
      ip_address INET
    )`,

    // Performance metrics table
    `CREATE TABLE IF NOT EXISTS performance_metrics (
      id SERIAL PRIMARY KEY,
      session_id VARCHAR(255) REFERENCES sessions(id),
      fps INTEGER,
      frame_time DECIMAL,
      memory_usage DECIMAL,
      triangles_in_view INTEGER,
      timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )`
  ]

  for (const query of queries) {
    try {
      await postgresPool.query(query)
    } catch (error) {
      logger.error('Failed to create table:', error)
      throw error
    }
  }

  logger.info('Database tables created/verified successfully')
}

// Main initialization function
export async function initializeDatabase() {
  try {
    await initializePostgreSQL()
    await initializeRedis() 
    await initializeNeo4j()
    logger.info('All databases initialized successfully')
  } catch (error) {
    logger.error('Database initialization failed:', error)
    throw error
  }
}

// Export database clients
export const getPostgresPool = () => postgresPool
export const getRedisClient = () => redisClient
export const getNeo4jDriver = () => neo4jDriver

// Graceful shutdown
export async function closeDatabaseConnections() {
  try {
    if (postgresPool) {
      await postgresPool.end()
      logger.info('PostgreSQL connection closed')
    }
    
    if (redisClient) {
      await redisClient.quit()
      logger.info('Redis connection closed')
    }
    
    if (neo4jDriver) {
      await neo4jDriver.close()
      logger.info('Neo4j connection closed')
    }
  } catch (error) {
    logger.error('Error closing database connections:', error)
  }
} 