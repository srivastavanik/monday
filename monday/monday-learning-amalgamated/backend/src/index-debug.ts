console.log('[DEBUG] Starting debug version...');

try {
  console.log('[DEBUG] Loading dotenv...');
  import('dotenv').then(dotenv => {
    import('path').then(path => {
      dotenv.config({ path: path.resolve(process.cwd(), '../.env') });
      console.log('[DEBUG] Environment loaded');
      console.log('[DEBUG] PERPLEXITY_API_KEY exists:', !!process.env.PERPLEXITY_API_KEY);
    });
  });
} catch (error) {
  console.error('[DEBUG] Error loading dotenv:', error);
}

try {
  console.log('[DEBUG] Loading express...');
  import('express').then(express => {
    console.log('[DEBUG] Express loaded');
    
    const app = express.default();
    const PORT = 3001;
    
    app.get('/health', (req, res) => {
      res.json({ status: 'debug-ok', timestamp: new Date().toISOString() });
    });
    
    app.listen(PORT, () => {
      console.log(`[DEBUG] Simple server running on port ${PORT}`);
    });
  });
} catch (error) {
  console.error('[DEBUG] Error loading express:', error);
}

console.log('[DEBUG] Debug script completed'); 