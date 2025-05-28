import type { Metadata } from 'next'
import './globals.css'

export const metadata: Metadata = {
  title: 'Monday - AI Learning Assistant',
  description: 'Voice-controlled AI assistant powered by Perplexity Sonar',
  generator: 'v0.dev',
  icons: {
    icon: '/favicon.svg',
  },
}

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  )
}
