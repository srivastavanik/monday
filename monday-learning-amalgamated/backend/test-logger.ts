console.log('Testing logger import...');

try {
  import('./src/utils/logger.js').then((module) => {
    console.log('Logger imported successfully:', typeof module.logger);
  }).catch((error) => {
    console.error('Logger import failed:', error);
  });
} catch (error) {
  console.error('Sync error:', error);
}

console.log('Test script completed'); 