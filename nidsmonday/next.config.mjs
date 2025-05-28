/** @type {import('next').NextConfig} */
const nextConfig = {
  eslint: {
    ignoreDuringBuilds: true,
  },
  typescript: {
    ignoreBuildErrors: true,
  },
  images: {
    unoptimized: true,
  },
  env: {
    NEXT_PUBLIC_ELEVENLABS_API_KEY: process.env.NEXT_PUBLIC_ELEVENLABS_API_KEY || 'sk_9454bca5ee475d45dc6e50bcb33bc3fb76f138e2191ff47d',
    NEXT_PUBLIC_ELEVENLABS_VOICE_ID: process.env.NEXT_PUBLIC_ELEVENLABS_VOICE_ID || '21m00Tcm4TlvDq8ikWAM',
    NEXT_PUBLIC_MONDAY_WS_URL: process.env.NEXT_PUBLIC_MONDAY_WS_URL || 'http://localhost:3001',
    NEXT_PUBLIC_MONDAY_API_URL: process.env.NEXT_PUBLIC_MONDAY_API_URL || 'http://localhost:3001/api'
  }
}

export default nextConfig
