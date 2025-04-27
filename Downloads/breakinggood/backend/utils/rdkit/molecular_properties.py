#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import base64
from io import BytesIO

def calculate_properties(smiles):
    """Calculate molecular properties for a given SMILES string"""
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}
        
        # Calculate basic properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        rings = Descriptors.RingCount(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        qed = QED.qed(mol)
        
        # Lipinski's Rule of Five
        lipinski_violations = Lipinski.NumRotatableBonds(mol) > 10
        lipinski_violations += mw > 500
        lipinski_violations += logp > 5
        lipinski_violations += hba > 10
        lipinski_violations += hbd > 5
        
        # Generate 2D depiction for visualization
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png_data = drawer.GetDrawingText()
        png_base64 = base64.b64encode(png_data).decode('utf-8')
        
        # Calculate fingerprints for similarity searching
        morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        morgan_fp_bits = [i for i in range(2048) if morgan_fp.GetBit(i)]
        
        # Synthetic accessibility score (lower is better)
        try:
            from rdkit.Chem import SynthesisHelpers
            synth_score = SynthesisHelpers.getScorefromSmiles(smiles)
        except:
            # If SynthesisHelpers not available, estimate based on complexity
            synth_score = 5.0 + (heavy_atoms / 100.0) + (rings / 10.0) * (rotatable_bonds / 10.0)
            if synth_score > 10:
                synth_score = 10.0
        
        # Estimate receptor binding affinity for ADHD targets (simplified model)
        # In a real implementation, this would use a more sophisticated ML model
        # trained on actual binding data for dopamine and norepinephrine transporters
        
        # Simple heuristic model based on known structure-activity relationships
        has_basic_nitrogen = 'N' in smiles and not 'NO' in smiles
        has_aromatic = '[nH]' in smiles or 'c1' in smiles
        has_phenyl = 'c1ccccc1' in smiles
        has_hydroxyl = 'O' in smiles # Simplified check for hydroxyl/ether
        has_halogen = any(x in smiles for x in ['F', 'Cl', 'Br', 'I'])

        dat_affinity = 0.0  # Dopamine transporter
        net_affinity = 0.0  # Norepinephrine transporter
        sert_affinity = 0.0 # Serotonin transporter
        d1_affinity = 0.0   # D1 receptor
        d2_affinity = 0.0   # D2 receptor

        if has_basic_nitrogen and has_aromatic:
            dat_affinity += 6.5
            net_affinity += 6.0
            sert_affinity += 5.0 # SERT often likes similar features but maybe less strongly
            d2_affinity += 6.0   # D2 often binds basic amines

        if has_phenyl:
            dat_affinity += 1.5
            net_affinity += 1.0
            sert_affinity += 1.8 # Phenyl groups can contribute to SERT binding
            d2_affinity += 1.0   # Phenyl interaction with D2

        if has_hydroxyl:
            net_affinity += 0.5  # Hydroxyls can interact with NET
            sert_affinity += 0.8 # And SERT
            d1_affinity += 1.5   # Catechol-like features are common for D1

        if has_halogen:
            sert_affinity += 1.0 # Halogens can sometimes enhance SERT binding
            d2_affinity += 0.5   # And D2

        if 200 < mw < 350:
            dat_affinity += 1.0
            net_affinity += 1.0
            sert_affinity += 0.8
            d1_affinity += 1.0
            d2_affinity += 1.2

        if 2 < logp < 4:
            dat_affinity += 1.0
            net_affinity += 1.5
            sert_affinity += 1.2
            d1_affinity += 0.5
            d2_affinity += 1.0

        # Randomize slightly to simulate model variance
        np.random.seed(int(sum(bytearray(smiles, 'utf-8'))))  # Deterministic based on SMILES
        dat_affinity += np.random.normal(0, 0.5)
        net_affinity += np.random.normal(0, 0.5)
        sert_affinity += np.random.normal(0, 0.6) # Slightly more variance
        d1_affinity += np.random.normal(0, 0.7)
        d2_affinity += np.random.normal(0, 0.5)

        # Clamp values to reasonable range (pKi values typically 4-10)
        dat_affinity = max(4.0, min(10.0, dat_affinity))
        net_affinity = max(4.0, min(10.0, net_affinity))
        sert_affinity = max(4.0, min(10.0, sert_affinity))
        d1_affinity = max(4.0, min(10.0, d1_affinity))
        d2_affinity = max(4.0, min(10.0, d2_affinity))

        # Create selectivity ratio (example: DAT/SERT)
        selectivity_dat_sert = dat_affinity / sert_affinity if sert_affinity > 0 else 1.0

        # Calculate additional CNS-relevant properties
        penetrates_bbb = (logp > 2.0 and logp < 5.0 and mw < 450 and tpsa < 90)
        cns_mpo_score = calculate_cns_mpo_score(mol)
        
        # Return calculated properties
        return {
            'smiles': smiles,
            'molecular_weight': round(mw, 2),
            'logp': round(logp, 2),
            'tpsa': round(tpsa, 2),
            'hba': hba,
            'hbd': hbd,
            'rotatable_bonds': rotatable_bonds,
            'rings': rings,
            'heavy_atoms': heavy_atoms,
            'qed': round(qed, 3),
            'lipinski_violations': lipinski_violations,
            'png_base64': png_base64,
            'synthetic_accessibility': round(synth_score, 2),
            'dat_affinity': round(dat_affinity, 2),
            'net_affinity': round(net_affinity, 2),
            'sert_affinity': round(sert_affinity, 2),
            'd1_affinity': round(d1_affinity, 2),
            'd2_affinity': round(d2_affinity, 2),
            'selectivity_dat_sert': round(selectivity_dat_sert, 2),
            'penetrates_bbb': penetrates_bbb,
            'cns_mpo_score': round(cns_mpo_score, 2),
            'morgan_fingerprint': morgan_fp_bits
        }
    except Exception as e:
        return {'error': str(e)}

def calculate_cns_mpo_score(mol):
    """Calculate CNS MPO score (0-6 scale) - higher is better for CNS drugs"""
    try:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        pka = 10  # Would need proper pKa calculator - this is placeholder
        
        # CNS MPO parameters and desirability functions
        # Based on ACS Chem Neurosci. 2010, 1, 435-449
        
        # Molecular weight (lower is better)
        d_mw = 1.0
        if mw > 360:
            d_mw = max(0, 1.0 - (mw - 360) / 40)
        
        # cLogP (2-5 is ideal)
        d_logp = 0.0
        if 2 <= logp <= 5:
            d_logp = 1.0
        elif logp < 2:
            d_logp = 0.25 * (logp + 2)
        else:
            d_logp = max(0, 1.0 - (logp - 5) / 2)
        
        # TPSA (40-90 is ideal)
        d_tpsa = 0.0
        if 40 <= tpsa <= 90:
            d_tpsa = 1.0
        elif tpsa < 40:
            d_tpsa = tpsa / 40
        else:
            d_tpsa = max(0, 1.0 - (tpsa - 90) / 50)
        
        # HBD (0-3 is ideal)
        d_hbd = 0.0
        if hbd <= 3:
            d_hbd = 1.0 - (hbd / 4)
        
        # Basic pKa (8-9.5 is ideal)
        d_pka = 0.0
        if 8 <= pka <= 9.5:
            d_pka = 1.0
        elif pka < 8:
            d_pka = pka / 8
        else:
            d_pka = max(0, 1.0 - (pka - 9.5) / 4.5)
        
        # Acidic pKa (placeholder - would need proper calculation)
        d_acidic = 1.0
        
        # Overall CNS MPO score
        cns_mpo = d_mw + d_logp + d_tpsa + d_hbd + d_pka + d_acidic
        return cns_mpo
        
    except Exception as e:
        print(f"Error calculating CNS MPO score: {str(e)}", file=sys.stderr)
        return 3.0  # Default middle value

def main():
    """Main function to process SMILES from command line arguments"""
    if len(sys.argv) < 2:
        print(json.dumps({'error': 'No SMILES provided'}))
        sys.exit(1)
    
    smiles = sys.argv[1]
    result = calculate_properties(smiles)
    print(json.dumps(result))

if __name__ == '__main__':
    main()
