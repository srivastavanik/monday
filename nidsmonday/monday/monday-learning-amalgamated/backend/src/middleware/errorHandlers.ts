// Custom error class for API errors
export class ApiError extends Error {
  public statusCode: number
  public isOperational: boolean

  constructor(message: string, statusCode: number = 500, isOperational: boolean = true) {
    super(message)
    this.statusCode = statusCode
    this.isOperational = isOperational
  }
}

// Error handler middleware
export const errorHandler = (error: any, req: any, res: any, next: any) => {
  let statusCode = 500
  let message = 'Internal Server Error'

  if (error instanceof ApiError) {
    statusCode = error.statusCode
    message = error.message
  }

  res.status(statusCode).json({
    success: false,
    message,
    statusCode,
    timestamp: new Date().toISOString()
  })
}

// 404 handler
export const notFoundHandler = (req: any, res: any) => {
  const message = `Route ${req.method} ${req.url} not found`
  
  res.status(404).json({
    success: false,
    message,
    statusCode: 404,
    timestamp: new Date().toISOString()
  })
} 