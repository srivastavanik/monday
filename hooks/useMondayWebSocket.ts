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
    if (!socketUrl || retryCountRef.current >= maxRetries) {
      if (retryCountRef.current >= maxRetries) {
        console.warn('Socket.IO: Maximum retry attempts reached.');
      }
      return;
    }

    const socket = io(socketUrl, {
      transports: ['websocket', 'polling'],
      timeout: 20000,
    });
    
    socketRef.current = socket;

    socket.on('connect', () => {
      console.log('Socket.IO: Connection established');
      setIsConnected(true);
      setError(null);
      setReadyState(1); // OPEN
      retryCountRef.current = 0;
      if (onOpen) onOpen(new Event('open'));
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

    socket.on('connect_error', (err) => {
      console.error('Socket.IO: Connection error', err);
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