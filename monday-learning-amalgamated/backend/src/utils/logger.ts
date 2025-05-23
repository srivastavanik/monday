import winston from 'winston'

const logLevel = process.env.LOG_LEVEL || 'info'
const isDevelopment = process.env.NODE_ENV === 'development'

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
      log += ` ${JSON.stringify(meta, null, 2)}`
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

// Add error logging helper
export const logError = (error: Error, context?: Record<string, any>) => {
  logger.error('Application Error', {
    message: error.message,
    stack: error.stack,
    ...context
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