import Image from "next/image"

export function PerplexityLogo({ className = "w-8 h-8" }: { className?: string }) {
  return (
    <Image 
      src="/Perplexity-AI-App-Icon-2023.png (1).svg" 
      alt="Perplexity Logo" 
      width={32} 
      height={32} 
      className={className} 
    />
  )
} 