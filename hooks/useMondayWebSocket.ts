"use client"

import { useState, useEffect, useCallback, useRef } from 'react';
import { io, Socket } from 'socket.io-client';

interface WebSocketHookOptions {
  socketUrl: string;
  onOpen?: (event: Event) => void;
  onMessage?: (event: MessageEvent) => void;
  onClose?: (event: CloseEvent) => void;
  onError?: (event: Event) => void;
  retryInterval?: number;
  maxRetries?: number;
}

// Singleton socket manager
class SocketManager {
  private static instance: SocketManager;
  private socket: Socket | null = null;
  private listeners: Map<string, Set<Function>> = new Map();
  private connectionPromise: Promise<Socket> | null = null;

  private constructor() {}

  static getInstance(): SocketManager {
    if (!SocketManager.instance) {
      SocketManager.instance = new SocketManager();
      
      // Register for cleanup on hot reload (development only)
      if (process.env.NODE_ENV === 'development' && typeof window !== 'undefined') {
        // Store instance reference for hot reload cleanup
        (window as any).__mondaySocketManager = SocketManager.instance;
        
        // Clean up previous instance if it exists
        const oldInstance = (window as any).__mondaySocketManagerPrevious;
        if (oldInstance && oldInstance !== SocketManager.instance) {
          console.log('Socket.IO: Cleaning up previous singleton instance');
          oldInstance.disconnect();
        }
        
        // Mark current as previous for next hot reload
        (window as any).__mondaySocketManagerPrevious = SocketManager.instance;
      }
    }
    return SocketManager.instance;
  }

  static resetInstance(): void {
    if (SocketManager.instance) {
      SocketManager.instance.disconnect();
      SocketManager.instance = null as any;
    }
  }

  async connect(url: string): Promise<Socket> {
    // If already connecting, wait for that
    if (this.connectionPromise) {
      console.log('Socket.IO: Waiting for existing connection...');
      return this.connectionPromise;
    }

    // If already connected, return existing socket
    if (this.socket?.connected) {
      console.log('Socket.IO: Already connected');
      return this.socket;
    }

    // Create new connection
    this.connectionPromise = new Promise((resolve, reject) => {
      console.log('Socket.IO: Creating singleton connection to:', url);

      // Clean up any existing socket
      if (this.socket) {
        console.log('Socket.IO: Cleaning up previous socket');
        this.socket.removeAllListeners();
        this.socket.close();
        this.socket = null;
      }

      const socket = io(url, {
        transports: ['websocket', 'polling'],
        timeout: 10000,
        reconnection: true,
        reconnectionAttempts: 10,
        reconnectionDelay: 1000,
        reconnectionDelayMax: 5000,
        autoConnect: true,
        forceNew: false,
        upgrade: true,
        rememberUpgrade: true,
      });

      this.socket = socket;

      // Handle connection success
      socket.once('connect', () => {
        console.log('Socket.IO: ✅ Singleton connected!', socket.id);
        this.connectionPromise = null;
        this.notifyListeners('connect', new Event('connect'));
        resolve(socket);
      });

      // Handle connection error
      socket.once('connect_error', (err: any) => {
        // Don't log xhr poll errors as they're often transient
        if (!err.message.includes('xhr poll error')) {
          console.error('Socket.IO: Connection error:', err.message);
        }
        this.connectionPromise = null;
        reject(err);
      });

      // Setup persistent event handlers
      socket.on('disconnect', (reason) => {
        console.log('Socket.IO: Disconnected:', reason);
        this.notifyListeners('disconnect', new CloseEvent('close', { reason }));
      });

      socket.on('connect_error', (err: any) => {
        // Only log serious errors, not transient polling errors
        if (err.type !== 'TransportError' && !err.message.includes('xhr poll error')) {
          console.error('Socket.IO: Error:', err.message);
        }
        this.notifyListeners('error', new Event('error'));
      });

      socket.on('reconnect', (attemptNumber) => {
        console.log('Socket.IO: Reconnected after', attemptNumber, 'attempts');
        this.notifyListeners('connect', new Event('connect'));
      });

      // Handle all other events
      socket.onAny((eventName, ...args) => {
        if (!['connect', 'disconnect', 'connect_error', 'reconnect'].includes(eventName)) {
          const messageEvent = new MessageEvent('message', {
            data: JSON.stringify({ type: eventName, ...args[0] })
          });
          this.notifyListeners('message', messageEvent);
        }
      });

      // Timeout for initial connection
      setTimeout(() => {
        if (!socket.connected) {
          console.error('Socket.IO: Initial connection timeout');
          this.connectionPromise = null;
          reject(new Error('Connection timeout'));
        }
      }, 15000);
    });

    return this.connectionPromise;
  }

  getSocket(): Socket | null {
    return this.socket;
  }

  isConnected(): boolean {
    return this.socket?.connected || false;
  }

  send(data: string): void {
    if (this.socket?.connected) {
      try {
        const parsedData = JSON.parse(data);
        this.socket.emit(parsedData.type, parsedData);
      } catch (e) {
        console.error('Socket.IO: Failed to send:', e);
      }
    } else {
      console.error('Socket.IO: Not connected');
    }
  }

  addListener(event: string, callback: Function): void {
    if (!this.listeners.has(event)) {
      this.listeners.set(event, new Set());
    }
    this.listeners.get(event)!.add(callback);
  }

  removeListener(event: string, callback: Function): void {
    this.listeners.get(event)?.delete(callback);
  }

  private notifyListeners(event: string, data: any): void {
    this.listeners.get(event)?.forEach(callback => {
      try {
        callback(data);
      } catch (e) {
        console.error('Socket.IO: Listener error:', e);
      }
    });
  }

  disconnect(): void {
    if (this.socket) {
      console.log('Socket.IO: Manually disconnecting');
      this.socket.removeAllListeners();
      this.socket.close();
      this.socket = null;
      this.connectionPromise = null;
    }
  }
}

const useMondayWebSocket = ({
  socketUrl,
  onOpen,
  onMessage,
  onClose,
  onError,
}: WebSocketHookOptions) => {
  const [isConnected, setIsConnected] = useState(false);
  const [lastMessage, setLastMessage] = useState<MessageEvent | null>(null);
  const [error, setError] = useState<Event | null>(null);
  const [readyState, setReadyState] = useState<number>(0);
  const managerRef = useRef<SocketManager>(SocketManager.getInstance());

  // Create stable callbacks
  const handleConnect = useCallback((event: Event) => {
    setIsConnected(true);
    setError(null);
    setReadyState(1);
    onOpen?.(event);
  }, [onOpen]);

  const handleDisconnect = useCallback((event: CloseEvent) => {
    setIsConnected(false);
    setReadyState(3);
    onClose?.(event);
  }, [onClose]);

  const handleError = useCallback((event: Event) => {
    setError(event);
    setReadyState(3);
    onError?.(event);
  }, [onError]);

  const handleMessage = useCallback((event: MessageEvent) => {
    setLastMessage(event);
    onMessage?.(event);
  }, [onMessage]);

  // Initialize connection
  useEffect(() => {
    const manager = managerRef.current;
    let mounted = true;

    const init = async () => {
      try {
        await manager.connect(socketUrl);
        if (mounted) {
          setIsConnected(manager.isConnected());
          setReadyState(manager.isConnected() ? 1 : 0);
        }
      } catch (err) {
        console.error('Socket.IO: Connection failed:', err);
        if (mounted) {
          setIsConnected(false);
          setError(new Event('error'));
          setReadyState(3);
        }
      }
    };

    // Add listeners
    manager.addListener('connect', handleConnect);
    manager.addListener('disconnect', handleDisconnect);
    manager.addListener('error', handleError);
    manager.addListener('message', handleMessage);

    // Initialize connection
    init();

    // Cleanup
    return () => {
      mounted = false;
      manager.removeListener('connect', handleConnect);
      manager.removeListener('disconnect', handleDisconnect);
      manager.removeListener('error', handleError);
      manager.removeListener('message', handleMessage);
      // Don't disconnect - let the singleton persist
    };
  }, [socketUrl, handleConnect, handleDisconnect, handleError, handleMessage]);

  // Monitor connection state
  useEffect(() => {
    const interval = setInterval(() => {
      const manager = managerRef.current;
      const connected = manager.isConnected();
      if (connected !== isConnected) {
        console.log('Socket.IO: State sync - updating to:', connected);
        setIsConnected(connected);
        setReadyState(connected ? 1 : 3);
      }
      
      // Send heartbeat if connected
      if (connected) {
        const socket = manager.getSocket();
        socket?.emit('heartbeat', { timestamp: Date.now() });
      }
    }, 10000);

    return () => clearInterval(interval);
  }, [isConnected]);

  const sendMessage = useCallback((data: string) => {
    managerRef.current.send(data);
  }, []);

  const reconnect = useCallback(async () => {
    const manager = managerRef.current;
    manager.disconnect();
    try {
      await manager.connect(socketUrl);
      setIsConnected(manager.isConnected());
      setReadyState(manager.isConnected() ? 1 : 0);
    } catch (err) {
      console.error('Socket.IO: Reconnection failed:', err);
      setIsConnected(false);
      setError(new Event('error'));
      setReadyState(3);
    }
  }, [socketUrl]);

  // Debug logging
  useEffect(() => {
    if (process.env.NODE_ENV === 'development') {
      const interval = setInterval(() => {
        const manager = managerRef.current;
        const socket = manager.getSocket();
        console.log('Socket.IO Debug:', {
          isConnected,
          managerConnected: manager.isConnected(),
          socketId: socket?.id,
          readyState,
        });
      }, 30000); // Every 30 seconds instead of 5

      return () => clearInterval(interval);
    }
  }, [isConnected, readyState]);

  return {
    isConnected,
    lastMessage,
    error,
    sendMessage,
    readyState,
    reconnect,
  };
};

export default useMondayWebSocket; 