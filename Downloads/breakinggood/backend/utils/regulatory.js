const fs = require('fs');
const path = require('path');
const axios = require('axios');
const logger = require('./logger');

// Regulatory timeline constants (in months)
const REGULATORY_TIMELINES = {
  preclinical: {
    min: 12,
    avg: 18,
    max: 30,
    description: 'Laboratory and animal testing to determine safety and efficacy before human trials'
  },
  phase1: {
    min: 6,
    avg: 12,
    max: 18,
    description: 'Small group testing (20-100 healthy volunteers) to determine safety and dosage'
  },
  phase2: {
    min: 12,
    avg: 24,
    max: 36,
    description: 'Medium group testing (100-500 volunteer patients) to determine efficacy and side effects'
  },
  phase3: {
    min: 24,
    avg: 36,
    max: 48,
    description: 'Large group testing (1,000-5,000 volunteer patients) to confirm efficacy and monitor adverse reactions'
  },
  nda_review: {
    min: 8,
    avg: 12,
    max: 16,
    description: 'FDA review of New Drug Application (NDA)'
  },
  manufacturing_setup: {
    min: 12,
    avg: 18,
    max: 24,
    description: 'Setting up manufacturing facilities and processes'
  },
  post_approval: {
    min: 36,
    avg: 48,
    max: 60,
    description: 'Post-approval monitoring and additional studies (Phase 4)'
  }
};

/**
 * Calculate regulatory timeline for a candidate molecule
 * @param {Object} moleculeProperties - Properties of the molecule
 * @param {Object} admetData - ADMET predictions for the molecule
 * @param {Object} productionData - Production feasibility data
 * @returns {Object} Regulatory timeline and analysis
 */
const calculateRegulatoryTimeline = (moleculeProperties, admetData, productionData) => {
  try {
    // Determine risk factors that may extend timeline
    const riskFactors = [];
    let riskScore = 0;
    
    // ADMET risk factors
    if (admetData) {
      // Check for toxicity concerns
      if (admetData.toxicity.herg_inhibition === 'high') {
        riskFactors.push('High hERG inhibition risk may require additional cardiac safety studies');
        riskScore += 0.15;
      }
      
      if (admetData.toxicity.hepatotoxicity === 'high') {
        riskFactors.push('High hepatotoxicity risk may require additional liver safety monitoring');
        riskScore += 0.15;
      }
      
      if (admetData.toxicity.cardiotoxicity_risk > 0.5) {
        riskFactors.push('Elevated cardiotoxicity risk will require thorough QT studies');
        riskScore += 0.1;
      }
      
      if (admetData.toxicity.ames_toxic) {
        riskFactors.push('Potential mutagenicity flagged by AMES prediction requires additional genotoxicity testing');
        riskScore += 0.2;
      }
      
      // Check for metabolism concerns
      if (admetData.metabolism.cyp_inhibition.CYP3A4 === 'yes') {
        riskFactors.push('CYP3A4 inhibition increases risk of drug-drug interactions');
        riskScore += 0.05;
      }
      
      if (admetData.metabolism.cyp_inhibition.CYP2D6 === 'yes') {
        riskFactors.push('CYP2D6 inhibition increases risk of drug-drug interactions');
        riskScore += 0.05;
      }
      
      // Check for addiction potential
      if (admetData.toxicity.addiction_potential === 'high') {
        riskFactors.push('High addiction potential will trigger DEA scheduling and additional abuse liability studies');
        riskScore += 0.25;
      } else if (admetData.toxicity.addiction_potential === 'medium') {
        riskFactors.push('Medium addiction potential may require abuse liability assessment');
        riskScore += 0.1;
      }
    }
    
    // Production feasibility risk factors
    if (productionData) {
      if (productionData.synthetic_complexity > 7) {
        riskFactors.push('Complex synthesis may create manufacturing challenges and delays');
        riskScore += 0.15;
      }
      
      if (productionData.scale_up_risk === 'high') {
        riskFactors.push('High scale-up risk may require additional process development time');
        riskScore += 0.1;
      }
      
      if (productionData.estimated_cost_per_kg > 10000) {
        riskFactors.push('High production cost may impact commercial viability');
        riskScore += 0.05;
      }
    }
    
    // Molecular property risk factors
    if (moleculeProperties) {
      if (moleculeProperties.lipinski_violations > 1) {
        riskFactors.push('Multiple Lipinski rule violations may impact drug development success');
        riskScore += 0.1;
      }
      
      if (moleculeProperties.synthetic_accessibility > 7) {
        riskFactors.push('Synthetic accessibility challenges may delay development');
        riskScore += 0.1;
      }
    }
    
    // Calculate timeline with risk adjustments
    const timeline = {};
    const startDate = new Date();
    let currentDate = new Date(startDate);
    
    // Preclinical phase
    const preclinical = { ...REGULATORY_TIMELINES.preclinical };
    preclinical.start_date = new Date(currentDate);
    // Adjust timeline based on risk factors
    const preclinical_adjustment = Math.round(preclinical.avg * riskScore);
    preclinical.adjusted_duration = preclinical.avg + preclinical_adjustment;
    currentDate.setMonth(currentDate.getMonth() + preclinical.adjusted_duration);
    preclinical.end_date = new Date(currentDate);
    timeline.preclinical = preclinical;
    
    // Phase 1
    const phase1 = { ...REGULATORY_TIMELINES.phase1 };
    phase1.start_date = new Date(currentDate);
    const phase1_adjustment = Math.round(phase1.avg * riskScore);
    phase1.adjusted_duration = phase1.avg + phase1_adjustment;
    currentDate.setMonth(currentDate.getMonth() + phase1.adjusted_duration);
    phase1.end_date = new Date(currentDate);
    timeline.phase1 = phase1;
    
    // Phase 2
    const phase2 = { ...REGULATORY_TIMELINES.phase2 };
    phase2.start_date = new Date(currentDate);
    const phase2_adjustment = Math.round(phase2.avg * riskScore);
    phase2.adjusted_duration = phase2.avg + phase2_adjustment;
    currentDate.setMonth(currentDate.getMonth() + phase2.adjusted_duration);
    phase2.end_date = new Date(currentDate);
    timeline.phase2 = phase2;
    
    // Phase 3
    const phase3 = { ...REGULATORY_TIMELINES.phase3 };
    phase3.start_date = new Date(currentDate);
    const phase3_adjustment = Math.round(phase3.avg * riskScore);
    phase3.adjusted_duration = phase3.avg + phase3_adjustment;
    currentDate.setMonth(currentDate.getMonth() + phase3.adjusted_duration);
    phase3.end_date = new Date(currentDate);
    timeline.phase3 = phase3;
    
    // Manufacturing can happen in parallel with Phase 3
    const manufacturing = { ...REGULATORY_TIMELINES.manufacturing_setup };
    manufacturing.start_date = new Date(phase3.start_date);
    const manufacturing_adjustment = Math.round(manufacturing.avg * (productionData ? productionData.scale_up_risk_factor : riskScore));
    manufacturing.adjusted_duration = manufacturing.avg + manufacturing_adjustment;
    const manufacturingEndDate = new Date(manufacturing.start_date);
    manufacturingEndDate.setMonth(manufacturingEndDate.getMonth() + manufacturing.adjusted_duration);
    manufacturing.end_date = manufacturingEndDate;
    timeline.manufacturing_setup = manufacturing;
    
    // NDA Review
    const nda_review = { ...REGULATORY_TIMELINES.nda_review };
    nda_review.start_date = new Date(currentDate);
    const nda_adjustment = Math.round(nda_review.avg * (riskScore / 2)); // Less impact on regulatory review time
    nda_review.adjusted_duration = nda_review.avg + nda_adjustment;
    currentDate.setMonth(currentDate.getMonth() + nda_review.adjusted_duration);
    nda_review.end_date = new Date(currentDate);
    timeline.nda_review = nda_review;
    
    // Post-approval monitoring
    const post_approval = { ...REGULATORY_TIMELINES.post_approval };
    post_approval.start_date = new Date(currentDate);
    post_approval.adjusted_duration = post_approval.avg; // Not adjusted by risk
    currentDate.setMonth(currentDate.getMonth() + post_approval.adjusted_duration);
    post_approval.end_date = new Date(currentDate);
    timeline.post_approval = post_approval;
    
    // Calculate total time to market (TTM)
    const ttm = {
      min: REGULATORY_TIMELINES.preclinical.min + 
           REGULATORY_TIMELINES.phase1.min + 
           REGULATORY_TIMELINES.phase2.min + 
           REGULATORY_TIMELINES.phase3.min + 
           REGULATORY_TIMELINES.nda_review.min,
      
      avg: REGULATORY_TIMELINES.preclinical.avg + 
           REGULATORY_TIMELINES.phase1.avg + 
           REGULATORY_TIMELINES.phase2.avg + 
           REGULATORY_TIMELINES.phase3.avg + 
           REGULATORY_TIMELINES.nda_review.avg,
      
      max: REGULATORY_TIMELINES.preclinical.max + 
           REGULATORY_TIMELINES.phase1.max + 
           REGULATORY_TIMELINES.phase2.max + 
           REGULATORY_TIMELINES.phase3.max + 
           REGULATORY_TIMELINES.nda_review.max,
      
      risk_adjusted: timeline.preclinical.adjusted_duration + 
                    timeline.phase1.adjusted_duration + 
                    timeline.phase2.adjusted_duration + 
                    timeline.phase3.adjusted_duration + 
                    timeline.nda_review.adjusted_duration
    };
    
    // Calculate approval probability
    let approvalProbability = 0.1; // Base 10% chance to account for preclinical success
    
    // Success rates for each phase (industry averages for CNS drugs)
    const phaseSuccessRates = {
      phase1: 0.6,  // 60% of Phase 1 trials succeed
      phase2: 0.35, // 35% of Phase 2 trials succeed
      phase3: 0.55, // 55% of Phase 3 trials succeed
      nda: 0.85     // 85% of NDAs are approved
    };
    
    // Adjust success rates based on risk factors
    const adjustedSuccessRates = {
      phase1: Math.max(0.05, phaseSuccessRates.phase1 * (1 - riskScore)),
      phase2: Math.max(0.05, phaseSuccessRates.phase2 * (1 - riskScore)),
      phase3: Math.max(0.05, phaseSuccessRates.phase3 * (1 - riskScore)),
      nda: Math.max(0.2, phaseSuccessRates.nda * (1 - (riskScore / 2)))
    };
    
    // Cumulative probability of approval
    approvalProbability = 0.1 * 
                         adjustedSuccessRates.phase1 * 
                         adjustedSuccessRates.phase2 * 
                         adjustedSuccessRates.phase3 * 
                         adjustedSuccessRates.nda;
    
    // Return the regulatory analysis
    return {
      timeline,
      risk_factors: riskFactors,
      risk_score: parseFloat(riskScore.toFixed(2)),
      time_to_market: {
        months: ttm,
        years: {
          min: (ttm.min / 12).toFixed(1),
          avg: (ttm.avg / 12).toFixed(1),
          max: (ttm.max / 12).toFixed(1),
          risk_adjusted: (ttm.risk_adjusted / 12).toFixed(1)
        }
      },
      approval_probability: parseFloat((approvalProbability * 100).toFixed(1)),
      phase_success_rates: adjustedSuccessRates,
      estimated_approval_date: timeline.nda_review.end_date.toISOString().split('T')[0],
      estimated_launch_date: (() => {
        const launchDate = new Date(timeline.nda_review.end_date);
        launchDate.setMonth(launchDate.getMonth() + 3); // Typically 3 months after approval
        return launchDate.toISOString().split('T')[0];
      })()
    };
  } catch (error) {
    logger.error(`Error calculating regulatory timeline: ${error.message}`);
    return {
      error: error.message,
      timeline: REGULATORY_TIMELINES
    };
  }
};

module.exports = {
  REGULATORY_TIMELINES,
  calculateRegulatoryTimeline
};
