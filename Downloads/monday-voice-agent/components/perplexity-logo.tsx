import Image from "next/image"

export function PerplexityLogo({ className = "w-8 h-8" }: { className?: string }) {
  return <Image src="/perplexity-logo.png" alt="Perplexity Logo" width={32} height={32} className={className} />
}
