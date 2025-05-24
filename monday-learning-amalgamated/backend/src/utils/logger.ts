import winston from 'winston'

const logLevel = process.env.LOG_LEVEL || 'info'
const isDevelopment = process.env.NODE_ENV === 'development'

// Safe JSON stringify to handle circular references
const safeJsonStringify = (obj: any, indent?: number): string => {
  const seen = new WeakSet()
  return JSON.stringify(obj, (key, value) => {
    if (typeof value === 'object' && value !== null) {
      if (seen.has(value)) {
        return '[Circular Reference]'
      }
      seen.add(value)
    }
    // Filter out sensitive or problematic properties
    if (key === 'password' || key === 'token' || key === 'secret') {
      return '[REDACTED]'
    }
    // Handle Error objects specifically
    if (value instanceof Error) {
      return {
        name: value.name,
        message: value.message,
        stack: value.stack
      }
    }
    // Handle request/response objects
    if (key === 'req' || key === 'res' || key === 'request' || key === 'response') {
      return '[HTTP Object]'
    }
    return value
  }, indent)
}

// Custom log format
const logFormat = winston.format.combine(
  winston.format.timestamp({
    format: 'YYYY-MM-DD HH:mm:ss'
  }),
  winston.format.errors({ stack: true }),
  winston.format.json(),
  winston.format.prettyPrint()
)

// Console format for development
const consoleFormat = winston.format.combine(
  winston.format.colorize(),
  winston.format.timestamp({
    format: 'HH:mm:ss'
  }),
  winston.format.printf(({ timestamp, level, message, ...meta }) => {
    let log = `${timestamp} [${level}]: ${message}`
    if (Object.keys(meta).length > 0) {
      try {
        log += ` ${safeJsonStringify(meta, 2)}`
      } catch (error: any) {
        log += ` [Error serializing metadata: ${error.message}]`
      }
    }
    return log
  })
)

// Create the logger
export const logger = winston.createLogger({
  level: logLevel,
  format: logFormat,
  defaultMeta: { service: 'monday-backend' },
  transports: [
    // Console transport
    new winston.transports.Console({
      format: isDevelopment ? consoleFormat : logFormat,
      silent: process.env.NODE_ENV === 'test'
    }),
    
    // File transports for production
    ...(isDevelopment ? [] : [
      new winston.transports.File({
        filename: 'logs/error.log',
        level: 'error',
        maxsize: 5242880, // 5MB
        maxFiles: 5
      }),
      new winston.transports.File({
        filename: 'logs/combined.log',
        maxsize: 5242880, // 5MB
        maxFiles: 5
      })
    ])
  ],
  
  // Handle exceptions and rejections
  exceptionHandlers: isDevelopment ? [] : [
    new winston.transports.File({
      filename: 'logs/exceptions.log'
    })
  ],
  rejectionHandlers: isDevelopment ? [] : [
    new winston.transports.File({
      filename: 'logs/rejections.log'
    })
  ]
})

// Add request logging helper
export const logRequest = (req: any, res: any, responseTime: number) => {
  logger.info('HTTP Request', {
    method: req.method,
    url: req.url,
    statusCode: res.statusCode,
    responseTime: `${responseTime}ms`,
    userAgent: req.get('User-Agent'),
    ip: req.ip,
    contentLength: res.get('Content-Length')
  })
}

// Add error logging helper with safe serialization
export const logError = (error: Error, context?: Record<string, any>) => {
  const errorInfo = {
    name: error.name,
    message: error.message,
    stack: error.stack
  }
  
  const safeContext = context ? JSON.parse(safeJsonStringify(context)) : {}
  
  logger.error('Application Error', {
    error: errorInfo,
    ...safeContext
  })
}

// Add performance logging helper
export const logPerformance = (operation: string, duration: number, metadata?: Record<string, any>) => {
  logger.info('Performance Metric', {
    operation,
    duration: `${duration}ms`,
    ...metadata
  })
} 