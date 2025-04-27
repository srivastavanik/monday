#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors, AllChem, Crippen
import numpy as np

def predict_admet(smiles):
    """Predict ADMET properties for a given SMILES string"""
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}
        
        # Calculate basic descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        qed = QED.qed(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        
        # Absorption prediction
        # Caco-2 permeability (simplified model)
        caco2 = 0.0
        if (logp > 0 and mw < 500 and rotatable_bonds < 10 and tpsa < 120 and hbd < 5 and hba < 10):
            caco2 = -0.008 * tpsa + 0.05 * logp + 0.002 * mw + 0.5
            # Add some noise for realism
            np.random.seed(int(sum(bytearray(smiles, 'utf-8'))))
            caco2 = max(0, min(2, caco2 + np.random.normal(0, 0.1)))
        
        # F (Oral Bioavailability) prediction
        # Using a simplified rule-based approach
        f_oral = 0.0
        if (mw < 500 and logp > 0 and logp < 5 and rotatable_bonds < 10 and tpsa < 140):
            f_oral = 0.7 - 0.006 * tpsa - 0.005 * rotatable_bonds + 0.01 * qed
            f_oral = max(0, min(1, f_oral + np.random.normal(0, 0.05)))
        
        # Distribution prediction
        # Blood-Brain Barrier penetration
        logbb = 0.0
        bbb_penetration = False
        if (logp > 0 and mw < 450 and tpsa < 90 and hbd < 3):
            logbb = 0.02 * mw + 0.1 * logp - 0.01 * tpsa + 0.5
            logbb = max(-3, min(1.5, logbb + np.random.normal(0, 0.2)))
            bbb_penetration = logbb > -0.5
        
        # Plasma Protein Binding
        ppb = 0.5
        if (logp > 0):
            ppb = 0.2 + 0.1 * logp
            ppb = max(0.1, min(0.99, ppb + np.random.normal(0, 0.05)))
        
        # Metabolism prediction
        # CYP450 interactions simplified
        cyp_inhibition = {}
        # CYP2D6 prediction
        cyp_inhibition['CYP2D6'] = 'yes' if (logp > 3 and 'N' in smiles and mw > 200) else 'no'
        # CYP3A4 prediction
        cyp_inhibition['CYP3A4'] = 'yes' if (logp > 3.5 and mw > 250 and 'O' in smiles) else 'no'
        # CYP2C9 prediction
        cyp_inhibition['CYP2C9'] = 'yes' if (logp > 4 and mw > 230) else 'no'
        
        # Elimination prediction
        # Half-life estimation (very simplified)
        half_life_category = 'medium'
        if (mw > 400 and logp > 5):
            half_life_category = 'long'
        elif (mw < 300 and logp < 3):
            half_life_category = 'short'
        
        # Toxicity prediction
        # hERG Inhibition
        herg_risk = 'low'
        if (logp > 3.7 and 'N' in smiles and mw > 250):
            herg_risk = 'high'
        elif (2.5 < logp < 3.7 and 'N' in smiles):
            herg_risk = 'medium'
        
        # Hepatotoxicity risk
        hepatotox_risk = 'low'
        if (logp > 3 and ('c1ccc(N)cc1' in smiles or 'c1ccc(O)cc1' in smiles)):
            hepatotox_risk = 'medium'
        elif (logp > 5 or ('c1ccc(Cl)cc1' in smiles and logp > 3)):
            hepatotox_risk = 'high'
        
        # AMES mutagenicity prediction
        ames_toxic = False
        if ('c1ccc(N=N)cc1' in smiles or 'c1ccc(NCO)cc1' in smiles):
            ames_toxic = True
        
        # Cardiotoxicity risk prediction
        cardiotox_risk = 0.0
        if (herg_risk == 'high'):
            cardiotox_risk = 0.8
        elif (herg_risk == 'medium'):
            cardiotox_risk = 0.4
        else:
            # Basic estimate from molecular properties
            cardiotox_risk = 0.1 + 0.05 * logp + (0.001 * mw) - 0.2 * qed
            cardiotox_risk = max(0, min(1, cardiotox_risk))
        
        # Addiction/abuse potential prediction
        # Very simplified approach - would need a more sophisticated model in real life
        addiction_potential = 'low'
        if (bbb_penetration and 'c1ccccc1' in smiles):
            addiction_potential = 'medium'
        if (bbb_penetration and 'N' in smiles and tpsa < 70 and mw < 350):
            addiction_potential = 'high'
        
        # Return predicted ADMET properties
        return {
            'smiles': smiles,
            'absorption': {
                'caco2_permeability': round(caco2, 2),
                'oral_bioavailability': round(f_oral * 100, 1)
            },
            'distribution': {
                'bbb_penetration': bbb_penetration,
                'logbb': round(logbb, 2),
                'plasma_protein_binding': round(ppb * 100, 1)
            },
            'metabolism': {
                'cyp_inhibition': cyp_inhibition
            },
            'elimination': {
                'half_life': half_life_category
            },
            'toxicity': {
                'herg_inhibition': herg_risk,
                'hepatotoxicity': hepatotox_risk,
                'cardiotoxicity_risk': round(cardiotox_risk, 2),
                'ames_toxic': ames_toxic,
                'addiction_potential': addiction_potential
            },
            'overall_admet_score': calculate_admet_score({
                'caco2': caco2,
                'f_oral': f_oral,
                'bbb': bbb_penetration,
                'ppb': ppb,
                'cyp': cyp_inhibition,
                'herg': herg_risk,
                'hepatotox': hepatotox_risk,
                'cardiotox': cardiotox_risk,
                'ames': ames_toxic,
                'addiction': addiction_potential
            })
        }
    except Exception as e:
        return {'error': str(e)}

def calculate_admet_score(admet_properties):
    """Calculate an overall ADMET score from individual properties"""
    score = 0.0
    max_score = 10.0
    
    # Absorption - up to 2.0 points
    score += admet_properties['caco2'] / 2.0  # 0-1 pts
    score += admet_properties['f_oral']  # 0-1 pts
    
    # Distribution - up to 2.0 points
    if not admet_properties['bbb']:  # For ADHD drugs, we don't want excessive BBB penetration
        score += 0.5
    else:  # But some penetration is still good
        score += 1.0
    score += (1 - admet_properties['ppb'])  # Lower protein binding is better
    
    # Metabolism - up to 2.0 points
    cyp_score = 0.0
    for enzyme, inhibits in admet_properties['cyp'].items():
        if inhibits == 'no':
            cyp_score += 1.0
    score += (cyp_score / len(admet_properties['cyp'])) * 2.0
    
    # Toxicity - up to 4.0 points
    if admet_properties['herg'] == 'low':
        score += 1.0
    elif admet_properties['herg'] == 'medium':
        score += 0.5
    
    if admet_properties['hepatotox'] == 'low':
        score += 1.0
    elif admet_properties['hepatotox'] == 'medium':
        score += 0.5
    
    score += (1 - admet_properties['cardiotox']) * 1.0
    
    if not admet_properties['ames']:
        score += 0.5
    
    if admet_properties['addiction'] == 'low':
        score += 0.5
    
    # Normalize to a 0-10 scale
    normalized_score = (score / max_score) * 10.0
    return round(normalized_score, 1)

def main():
    """Main function to process SMILES from command line arguments"""
    if len(sys.argv) < 2:
        print(json.dumps({'error': 'No SMILES provided'}))
        sys.exit(1)
    
    smiles = sys.argv[1]
    result = predict_admet(smiles)
    print(json.dumps(result))

if __name__ == '__main__':
    main()
