console.log('[TEST] Starting minimal server test...');

import express from 'express';
import { createServer } from 'http';
import { Server } from 'socket.io';
import dotenv from 'dotenv';
import path from 'path';

// Load environment variables
dotenv.config({ path: path.resolve(process.cwd(), '../.env') });
console.log('[TEST] Environment loaded, PORT:', process.env.PORT);

const app = express();
const server = createServer(app);
const io = new Server(server, {
  cors: {
    origin: "http://localhost:3000",
    methods: ["GET", "POST"],
    credentials: true
  }
});

const PORT = process.env.PORT || 3001;

app.get('/health', (req, res) => {
  res.json({ status: 'ok', timestamp: new Date().toISOString() });
});

io.on('connection', (socket) => {
  console.log(`[TEST] Client connected: ${socket.id}`);
  socket.on('disconnect', () => {
    console.log(`[TEST] Client disconnected: ${socket.id}`);
  });
});

server.listen(PORT, () => {
  console.log(`[TEST] Minimal server running on port ${PORT}`);
});

console.log('[TEST] Setup completed'); 