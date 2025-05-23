export const authenticateSocket = (socket: any, next: any) => {
  // For now, allow all connections
  // In production, implement proper JWT validation
  next()
} 