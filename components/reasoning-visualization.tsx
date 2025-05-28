"use client"

import type { ReasoningStep } from "../types/monday"
import { Brain } from "lucide-react"

interface ReasoningVisualizationProps {
  steps: ReasoningStep[]
  currentStep: number
}

export function ReasoningVisualization({ steps, currentStep }: ReasoningVisualizationProps) {
  return (
    <div className="space-y-4">
      <div className="flex items-center gap-2 mb-6">
        <Brain className="w-5 h-5 text-blue-400" />
        <span className="text-blue-300 font-semibold">REASONING PROCESS</span>
      </div>

      <div className="relative">
        {/* Connecting line */}
        <div className="absolute left-6 top-8 bottom-0 w-0.5 bg-blue-500/30"></div>

        {steps.map((step, index) => (
          <div key={step.id} className="relative flex gap-4 mb-6">
            {/* Step indicator */}
            <div
              className={`relative z-10 w-12 h-12 rounded-full border-2 flex items-center justify-center text-sm font-bold transition-all duration-500 ${
                index <= currentStep
                  ? "bg-blue-500 border-blue-400 text-white"
                  : "bg-slate-800 border-slate-600 text-slate-400"
              }`}
            >
              {index + 1}
            </div>

            {/* Step content */}
            <div
              className={`flex-1 transition-all duration-500 ${index <= currentStep ? "opacity-100" : "opacity-40"}`}
            >
              <h4 className="text-white font-semibold mb-2">{step.step}</h4>
              <p className="text-slate-300 text-sm leading-relaxed mb-2">{step.reasoning}</p>

              {/* Confidence indicator */}
              <div className="flex items-center gap-2">
                <span className="text-xs text-slate-400">Confidence:</span>
                <div className="w-20 h-2 bg-slate-700 rounded-full overflow-hidden">
                  <div
                    className="h-full bg-gradient-to-r from-blue-500 to-green-400 transition-all duration-1000"
                    style={{ width: `${step.confidence * 100}%` }}
                  ></div>
                </div>
                <span className="text-xs text-slate-400">{Math.round(step.confidence * 100)}%</span>
              </div>
            </div>
          </div>
        ))}
      </div>
    </div>
  )
}
