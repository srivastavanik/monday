const fs = require('fs');
const path = require('path');
const axios = require('axios');
const logger = require('./logger');

/**
 * Calculate production feasibility for a candidate molecule
 * @param {Object} moleculeProperties - Properties of the molecule
 * @param {Object} admetData - ADMET predictions for the molecule
 * @returns {Object} Production feasibility analysis
 */
const calculateProductionFeasibility = (moleculeProperties, admetData) => {
  try {
    if (!moleculeProperties) {
      return { error: 'Molecule properties are required' };
    }
    
    // Extract relevant properties
    const syntheticAccessibility = moleculeProperties.synthetic_accessibility || 5;
    const mw = moleculeProperties.molecular_weight || 300;
    const logp = moleculeProperties.logp || 2;
    const tpsa = moleculeProperties.tpsa || 60;
    const rings = moleculeProperties.rings || 2;
    const rotatableBonds = moleculeProperties.rotatable_bonds || 5;
    
    // Calculate synthetic complexity (1-10 scale, higher is more complex)
    const syntheticComplexity = syntheticAccessibility;
    
    // Estimate number of synthetic steps
    const syntheticSteps = Math.round(2 + (syntheticComplexity / 2));
    
    // Calculate scale-up risk (0-1 scale)
    let scaleUpRiskFactor = 0.0;
    
    // Complex molecules are harder to scale
    scaleUpRiskFactor += syntheticComplexity / 20; // Contributes up to 0.5
    
    // More rings increase complexity
    scaleUpRiskFactor += (rings > 3) ? (rings - 3) * 0.05 : 0;
    
    // Very high or low LogP can create solubility issues in manufacturing
    if (logp < 0 || logp > 5) {
      scaleUpRiskFactor += 0.1;
    }
    
    // Molecules with high MW are often harder to purify at scale
    if (mw > 500) {
      scaleUpRiskFactor += (mw - 500) / 1000; // Adds 0.1 per 100 MW over 500
    }
    
    // Highly flexible molecules can cause crystallization issues
    if (rotatableBonds > 8) {
      scaleUpRiskFactor += (rotatableBonds - 8) * 0.02;
    }
    
    // Cap the risk factor at 1.0
    scaleUpRiskFactor = Math.min(1.0, scaleUpRiskFactor);
    
    // Convert risk factor to category
    let scaleUpRisk = 'low';
    if (scaleUpRiskFactor > 0.7) {
      scaleUpRisk = 'high';
    } else if (scaleUpRiskFactor > 0.4) {
      scaleUpRisk = 'medium';
    }
    
    // Estimate yield
    // Base yield per step (higher synthetic complexity = lower yield)
    const baseYieldPerStep = Math.max(50, 90 - (syntheticComplexity * 5));
    
    // Overall yield calculation (compound yields across steps)
    const overallYield = Math.pow(baseYieldPerStep / 100, syntheticSteps) * 100;
    
    // Estimate production cost
    // Base cost per step (USD)
    const baseCostPerStep = 500 + (syntheticComplexity * 200);
    
    // Factor in scale (cost per kg decreases with scale)
    const smallScaleCost = baseCostPerStep * syntheticSteps * (100 / overallYield);
    const mediumScaleCost = smallScaleCost * 0.7; // 30% reduction at medium scale
    const largeScaleCost = smallScaleCost * 0.5;  // 50% reduction at large scale
    
    // Risk-adjusted costs (uncertainty increases with risk)
    const riskPremium = 1 + scaleUpRiskFactor;
    
    // Calculate failure risk
    const riskOfFailure = {
      process_development: Math.min(0.95, scaleUpRiskFactor * 1.5),
      scale_up: Math.min(0.95, scaleUpRiskFactor * 1.2),
      regulatory_compliance: Math.min(0.8, scaleUpRiskFactor)
    };
    
    // Overall risk of failure (probability of at least one failure)
    const overallRiskOfFailure = 1 - (
      (1 - riskOfFailure.process_development) * 
      (1 - riskOfFailure.scale_up) * 
      (1 - riskOfFailure.regulatory_compliance)
    );
    
    // Analyze impurity profile and process reliability
    let impurityProfile = 'favorable';
    let processReliability = 'high';
    
    if (syntheticComplexity > 7) {
      impurityProfile = 'challenging';
      processReliability = 'low';
    } else if (syntheticComplexity > 5) {
      impurityProfile = 'moderate';
      processReliability = 'medium';
    }
    
    // Environmental impact
    const environmentalImpact = {
      solvent_usage: (syntheticSteps * 10) + (syntheticComplexity * 5), // Liters per kg
      e_factor: Math.min(100, syntheticSteps * syntheticComplexity), // kg waste per kg product
      energy_consumption: syntheticSteps * 5 // kWh per kg
    };
    
    // Estimate timeline (in months)
    const timeline = {
      process_development: Math.round(3 + (syntheticComplexity / 2)),
      pilot_scale: Math.round(2 + scaleUpRiskFactor * 6),
      commercial_scale: Math.round(3 + scaleUpRiskFactor * 9)
    };
    
    timeline.total = timeline.process_development + timeline.pilot_scale + timeline.commercial_scale;
    
    // Return the feasibility analysis
    return {
      synthetic_complexity: parseFloat(syntheticComplexity.toFixed(1)),
      synthetic_steps: syntheticSteps,
      scale_up_risk: scaleUpRisk,
      scale_up_risk_factor: parseFloat(scaleUpRiskFactor.toFixed(2)),
      estimated_yield: {
        per_step: parseFloat(baseYieldPerStep.toFixed(1)),
        overall: parseFloat(overallYield.toFixed(1))
      },
      estimated_cost_per_kg: {
        small_scale: Math.round(smallScaleCost * riskPremium),
        medium_scale: Math.round(mediumScaleCost * riskPremium),
        large_scale: Math.round(largeScaleCost * riskPremium)
      },
      risk_of_failure: {
        process_development: parseFloat(riskOfFailure.process_development.toFixed(2)),
        scale_up: parseFloat(riskOfFailure.scale_up.toFixed(2)),
        regulatory_compliance: parseFloat(riskOfFailure.regulatory_compliance.toFixed(2)),
        overall: parseFloat(overallRiskOfFailure.toFixed(2))
      },
      impurity_profile: impurityProfile,
      process_reliability: processReliability,
      environmental_impact: environmentalImpact,
      timeline_months: timeline
    };
  } catch (error) {
    logger.error(`Error calculating production feasibility: ${error.message}`);
    return {
      error: error.message
    };
  }
};

module.exports = {
  calculateProductionFeasibility
};
