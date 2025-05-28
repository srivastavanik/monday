console.log('TEST: Starting incremental import test');

import express from 'express';
console.log('TEST: Express import successful');

import { createServer } from 'http';
console.log('TEST: HTTP import successful');

import dotenv from 'dotenv';
import path from 'path';
dotenv.config({ path: path.resolve(process.cwd(), '../.env') });
console.log('TEST: Environment variables loaded from parent directory');
console.log('TEST: PORT from env:', process.env.PORT);
console.log('TEST: DATABASE_URL from env:', process.env.DATABASE_URL ? 'Present' : 'Missing');

console.log('TEST: Testing logger import...');
try {
  const { logger } = await import('./utils/logger.js');
  console.log('TEST: Logger import successful');
} catch (error) {
  console.log('TEST: Logger import FAILED:', (error as Error).message);
  process.exit(1);
}

console.log('TEST: Testing middleware imports...');
try {
  const { errorHandler, notFoundHandler } = await import('./middleware/errorHandlers.js');
  console.log('TEST: Middleware errorHandlers import successful');
} catch (error) {
  console.log('TEST: Middleware errorHandlers import FAILED:', (error as Error).message);
  process.exit(1);
}

try {
  const { authenticateSocket } = await import('./middleware/socketAuth.js');
  console.log('TEST: Middleware socketAuth import successful');
} catch (error) {
  console.log('TEST: Middleware socketAuth import FAILED:', (error as Error).message);
  process.exit(1);
}

console.log('TEST: Testing database import...');
try {
  const { initializeDatabase } = await import('./database/index.js');
  console.log('TEST: Database import successful');
} catch (error) {
  console.log('TEST: Database import FAILED:', (error as Error).message);
  process.exit(1);
}

console.log('TEST: Testing Redis service import...');
try {
  const { initializeRedis } = await import('./services/redis.js');
  console.log('TEST: Redis service import successful');
} catch (error) {
  console.log('TEST: Redis service import FAILED:', (error as Error).message);
  process.exit(1);
}

console.log('TEST: Testing route imports...');
try {
  const queryRoutes = await import('./routes/query.js');
  console.log('TEST: Query routes import successful');
} catch (error) {
  console.log('TEST: Query routes import FAILED:', (error as Error).message);
  process.exit(1);
}

try {
  const sessionRoutes = await import('./routes/session.js');
  console.log('TEST: Session routes import successful');
} catch (error) {
  console.log('TEST: Session routes import FAILED:', (error as Error).message);
  process.exit(1);
}

try {
  const analyticsRoutes = await import('./routes/analytics.js');
  console.log('TEST: Analytics routes import successful');
} catch (error) {
  console.log('TEST: Analytics routes import FAILED:', (error as Error).message);
  process.exit(1);
}

console.log('TEST: Testing socket handler imports...');
try {
  const { handleVoiceQuery } = await import('./sockets/voiceHandler.js');
  console.log('TEST: Voice handler import successful');
} catch (error) {
  console.log('TEST: Voice handler import FAILED:', (error as Error).message);
  process.exit(1);
}

try {
  const { handleSpatialCommands } = await import('./sockets/spatialHandler.js');
  console.log('TEST: Spatial handler import successful');
} catch (error) {
  console.log('TEST: Spatial handler import FAILED:', (error as Error).message);
  process.exit(1);
}

try {
  const { handleSessionEvents } = await import('./sockets/sessionHandler.js');
  console.log('TEST: Session handler import successful');
} catch (error) {
  console.log('TEST: Session handler import FAILED:', (error as Error).message);
  process.exit(1);
}

console.log('TEST: All imports successful!');
process.exit(0); 