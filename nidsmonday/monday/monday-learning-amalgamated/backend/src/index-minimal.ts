console.log('MINIMAL: Starting minimal index');

import express from 'express';
import { createServer } from 'http';
import dotenv from 'dotenv';
import path from 'path';

dotenv.config({ path: path.resolve(process.cwd(), '../.env') });
console.log('MINIMAL: Environment loaded, PORT:', process.env.PORT);

const app = express();
const server = createServer(app);
const PORT = process.env.PORT || 3001;

app.get('/health', (req, res) => {
  res.json({ status: 'ok', message: 'Minimal backend running' });
});

server.listen(PORT, () => {
  console.log(`MINIMAL: Server running on port ${PORT}`);
  console.log('MINIMAL: Try visiting http://localhost:3001/health');
});

console.log('MINIMAL: Setup complete'); 