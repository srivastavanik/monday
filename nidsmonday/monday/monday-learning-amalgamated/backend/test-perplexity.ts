console.log('Testing perplexity import...');

// Load environment variables first
import dotenv from 'dotenv';
import path from 'path';

dotenv.config({ path: path.resolve(process.cwd(), '../.env') });
console.log('PERPLEXITY_API_KEY exists:', !!process.env.PERPLEXITY_API_KEY);

try {
  import('./src/services/perplexity.js').then((module) => {
    console.log('Perplexity imported successfully:', typeof module.perplexityService);
  }).catch((error) => {
    console.error('Perplexity import failed:', error);
  });
} catch (error) {
  console.error('Sync error:', error);
}

console.log('Test script completed'); 