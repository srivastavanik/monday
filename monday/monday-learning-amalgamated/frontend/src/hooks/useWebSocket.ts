import { useState, useEffect, useRef } from 'react'
import { io, Socket } from 'socket.io-client'

interface UseWebSocketReturn {
  socket: Socket | null
  isConnected: boolean
  error: string | null
  connect: () => void
  disconnect: () => void
}

export const useWebSocketConnection = (): UseWebSocketReturn => {
  const [socket, setSocket] = useState<Socket | null>(null)
  const [isConnected, setIsConnected] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const reconnectTimeoutRef = useRef<NodeJS.Timeout | null>(null)

  const connect = () => {
    if (socket?.connected) return

    const newSocket = io('ws://localhost:3001', {
      transports: ['websocket'],
      autoConnect: true,
      reconnection: true,
      reconnectionAttempts: 5,
      reconnectionDelay: 1000,
    })

    newSocket.on('connect', () => {
      console.log('WebSocket connected')
      setIsConnected(true)
      setError(null)
    })

    newSocket.on('disconnect', () => {
      console.log('WebSocket disconnected')
      setIsConnected(false)
    })

    newSocket.on('connect_error', (err) => {
      console.error('WebSocket connection error:', err)
      setError('Connection failed')
      setIsConnected(false)
    })

    setSocket(newSocket)
  }

  const disconnect = () => {
    if (socket) {
      socket.disconnect()
      setSocket(null)
      setIsConnected(false)
    }
  }

  useEffect(() => {
    connect()

    return () => {
      if (reconnectTimeoutRef.current) {
        clearTimeout(reconnectTimeoutRef.current)
      }
      disconnect()
    }
  }, [])

  return {
    socket,
    isConnected,
    error,
    connect,
    disconnect
  }
} 