import express from 'express'
import { createServer } from 'http'
import { Server } from 'socket.io'
import cors from 'cors'
import helmet from 'helmet'
import compression from 'compression'
import rateLimit from 'express-rate-limit'
import dotenv from 'dotenv'

import { logger } from './utils/logger.js'
import { errorHandler, notFoundHandler } from './middleware/errorHandlers.js'
import { authenticateSocket } from './middleware/socketAuth.js'
import { initializeDatabase } from './database/index.js'
import { initializeRedis } from './services/redis.js'

// Import routes
import queryRoutes from './routes/query.js'
import sessionRoutes from './routes/session.js'
import analyticsRoutes from './routes/analytics.js'

// Import socket handlers
import { handleVoiceQuery } from './sockets/voiceHandler.js'
import { handleSpatialCommands } from './sockets/spatialHandler.js'
import { handleSessionEvents } from './sockets/sessionHandler.js'

// Load environment variables
dotenv.config()

const app = express()
const server = createServer(app)
const io = new Server(server, {
  cors: {
    origin: process.env.FRONTEND_URL || "http://localhost:3000",
    methods: ["GET", "POST"],
    credentials: true
  },
  transports: ['websocket', 'polling']
})

const PORT = process.env.PORT || 3001

// Security middleware
app.use(helmet({
  crossOriginEmbedderPolicy: false, // Needed for WebXR
  contentSecurityPolicy: {
    directives: {
      defaultSrc: ["'self'"],
      scriptSrc: ["'self'", "'unsafe-inline'", "'unsafe-eval'"], // WebXR needs unsafe-eval
      styleSrc: ["'self'", "'unsafe-inline'"],
      imgSrc: ["'self'", "data:", "https:"],
      connectSrc: ["'self'", "wss:", "ws:"],
      fontSrc: ["'self'", "https:"],
      objectSrc: ["'none'"],
      mediaSrc: ["'self'"],
      frameSrc: ["'none'"],
    },
  },
}))

// CORS configuration
app.use(cors({
  origin: process.env.FRONTEND_URL || "http://localhost:3000",
  credentials: true,
  methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
  allowedHeaders: ['Content-Type', 'Authorization', 'X-Requested-With']
}))

// Rate limiting
const limiter = rateLimit({
  windowMs: 15 * 60 * 1000, // 15 minutes
  max: 100, // Limit each IP to 100 requests per windowMs
  message: 'Too many requests from this IP, please try again later.',
  standardHeaders: true,
  legacyHeaders: false,
})
app.use('/api/', limiter)

// Basic middleware
app.use(compression())
app.use(express.json({ limit: '10mb' }))
app.use(express.urlencoded({ extended: true, limit: '10mb' }))

// Logging middleware
app.use((req, res, next) => {
  logger.info(`${req.method} ${req.url}`, {
    ip: req.ip,
    userAgent: req.get('User-Agent'),
    timestamp: new Date().toISOString()
  })
  next()
})

// Health check endpoint
app.get('/health', (req, res) => {
  res.status(200).json({
    status: 'healthy',
    timestamp: new Date().toISOString(),
    uptime: process.uptime(),
    memory: process.memoryUsage(),
    version: process.env.npm_package_version || '1.0.0'
  })
})

// API routes
app.use('/api/query', queryRoutes)
app.use('/api/session', sessionRoutes)
app.use('/api/analytics', analyticsRoutes)

// Socket.IO authentication middleware
io.use(authenticateSocket)

// Socket.IO connection handling
io.on('connection', (socket) => {
  logger.info(`Client connected: ${socket.id}`)
  
  // Register socket event handlers
  handleVoiceQuery(socket, io)
  handleSpatialCommands(socket, io)
  handleSessionEvents(socket, io)
  
  socket.on('disconnect', (reason) => {
    logger.info(`Client disconnected: ${socket.id}, reason: ${reason}`)
  })
  
  socket.on('error', (error) => {
    logger.error(`Socket error for ${socket.id}:`, error)
  })
})

// Error handling middleware (must be last)
app.use(notFoundHandler)
app.use(errorHandler)

// Initialize services and start server
async function startServer() {
  try {
    // Initialize database connections
    await initializeDatabase()
    logger.info('Database initialized successfully')
    
    // Initialize Redis
    await initializeRedis()
    logger.info('Redis initialized successfully')
    
    // Start the server
    server.listen(PORT, () => {
      logger.info(`Monday backend server running on port ${PORT}`)
      logger.info(`Environment: ${process.env.NODE_ENV || 'development'}`)
      logger.info(`Frontend URL: ${process.env.FRONTEND_URL || 'http://localhost:3000'}`)
    })
    
  } catch (error) {
    logger.error('Failed to start server:', error)
    process.exit(1)
  }
}

// Graceful shutdown handling
process.on('SIGTERM', async () => {
  logger.info('SIGTERM received, shutting down gracefully')
  server.close(() => {
    logger.info('Server closed')
    process.exit(0)
  })
})

process.on('SIGINT', async () => {
  logger.info('SIGINT received, shutting down gracefully')
  server.close(() => {
    logger.info('Server closed')
    process.exit(0)
  })
})

// Handle uncaught exceptions
process.on('uncaughtException', (error) => {
  logger.error('Uncaught Exception:', error)
  process.exit(1)
})

process.on('unhandledRejection', (reason, promise) => {
  logger.error('Unhandled Rejection at:', promise, 'reason:', reason)
  process.exit(1)
})

// Start the server
startServer() 