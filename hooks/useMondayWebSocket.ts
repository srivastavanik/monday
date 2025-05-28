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

interface WebSocketState {
  isConnected: boolean;
  lastMessage: MessageEvent | null;
  error: Event | null;
  readyState: number;
}

const useMondayWebSocket = ({
  socketUrl,
  onOpen,
  onMessage,
  onClose,
  onError,
  retryInterval = 5000,
  maxRetries = 10,
}: WebSocketHookOptions) => {
  const socketRef = useRef<Socket | null>(null);
  const [isConnected, setIsConnected] = useState(false);
  const [lastMessage, setLastMessage] = useState<MessageEvent | null>(null);
  const [error, setError] = useState<Event | null>(null);
  const [readyState, setReadyState] = useState<number>(0);
  const retryCountRef = useRef(0);

  const connect = useCallback(() => {
    console.log('Socket.IO: Attempting to connect to:', socketUrl);
    if (!socketUrl || retryCountRef.current >= maxRetries) {
      if (retryCountRef.current >= maxRetries) {
        console.warn('Socket.IO: Maximum retry attempts reached.');
      }
      return;
    }

    console.log('Socket.IO: Creating socket connection...');
    const socket = io(socketUrl, {
      transports: ['websocket', 'polling'],
      timeout: 20000,
      reconnection: true,
      reconnectionAttempts: 5,
      reconnectionDelay: 1000,
    });
    
    socketRef.current = socket;
    console.log('Socket.IO: Socket instance created');

    socket.on('connect', () => {
      console.log('Socket.IO: ✅ Connection event fired!');
      console.log('Socket.IO: Socket ID:', socket.id);
      console.log('Socket.IO: Connected status:', socket.connected);
      setIsConnected(true);
      console.log('Socket.IO: isConnected state set to TRUE');
      setError(null);
      setReadyState(1); // OPEN
      retryCountRef.current = 0;
      if (onOpen) {
        console.log('Socket.IO: Calling onOpen callback');
        onOpen(new Event('open'));
      }
    });

    socket.on('disconnect', (reason) => {
      console.log(`Socket.IO: Connection closed (Reason: ${reason})`);
      setIsConnected(false);
      setReadyState(3); // CLOSED
      if (onClose) onClose(new CloseEvent('close', { reason }));

      // Attempt to reconnect if not a clean close and under max retries
      if (reason !== 'io client disconnect' && retryCountRef.current < maxRetries) {
        retryCountRef.current++;
        console.log(`Socket.IO: Attempting to reconnect... (Attempt ${retryCountRef.current}/${maxRetries})`);
        setTimeout(connect, retryInterval * Math.pow(2, Math.min(retryCountRef.current - 1, 4)));
      }
    });

    socket.on('connect_error', (err: any) => {
      console.error('Socket.IO: Connection error', err);
      console.error('Socket.IO: Error type:', err.type);
      console.error('Socket.IO: Error message:', err.message);
      console.error('Socket.IO: Socket URL:', socketUrl);
      console.error('Socket.IO: Error details:', err);
      
      const errorEvent = new Event('error');
      setError(errorEvent);
      setReadyState(3); // CLOSED
      if (onError) onError(errorEvent);
    });

    // Listen for all events and convert to MessageEvent format
    socket.onAny((eventName, ...args) => {
      if (eventName !== 'connect' && eventName !== 'disconnect' && eventName !== 'connect_error') {
        const messageEvent = new MessageEvent('message', {
          data: JSON.stringify({ type: eventName, ...args[0] })
        });
        setLastMessage(messageEvent);
        if (onMessage) onMessage(messageEvent);
      }
    });

  }, [socketUrl, onOpen, onMessage, onClose, onError, retryInterval, maxRetries]);

  useEffect(() => {
    connect();

    return () => {
      if (socketRef.current) {
        console.log('Socket.IO: Cleaning up connection');
        socketRef.current.disconnect();
        socketRef.current = null;
      }
    };
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [socketUrl]);

  const sendMessage = useCallback((data: string) => {
    if (socketRef.current && socketRef.current.connected) {
      try {
        const parsedData = JSON.parse(data);
        socketRef.current.emit(parsedData.type, parsedData);
      } catch (e) {
        console.error('Socket.IO: Failed to parse message data', e);
      }
    } else {
      console.error('Socket.IO: Connection not open. Cannot send message.');
    }
  }, []);

  // Debug connection state every 2 seconds
  useEffect(() => {
    const interval = setInterval(() => {
      console.log('Socket.IO State Check:', {
        isConnected,
        socketExists: !!socketRef.current,
        socketConnected: socketRef.current?.connected,
        readyState,
        error: error ? 'Error present' : 'No error'
      });
    }, 2000);
    
    return () => clearInterval(interval);
  }, [isConnected, readyState, error]);

  return {
    isConnected,
    lastMessage,
    error,
    sendMessage,
    readyState,
    socketRef,
  };
};

export default useMondayWebSocket; 