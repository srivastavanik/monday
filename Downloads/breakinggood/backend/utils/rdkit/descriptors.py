#!/usr/bin/env python3
"""
Molecular descriptor calculator using RDKit.
This script calculates various molecular properties for drug discovery.
"""

import sys
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def calculate_descriptors(smiles):
    """Calculate molecular descriptors for a given SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": f"Invalid SMILES: {smiles}"}
        
        # Try to generate 3D coordinates for the molecule
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            mol = Chem.RemoveHs(mol)
            has_3d = True
        except:
            has_3d = False
        
        # Calculate basic properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        rings = Descriptors.RingCount(mol)
        aromatic_rings = sum(1 for ring in mol.GetSSSR() if ring.GetNumAtoms() >= 5 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
        heavy_atoms = mol.GetNumHeavyAtoms()
        
        # Lipinski's Rule of 5
        ro5_violations = 0
        if mw > 500: ro5_violations += 1
        if logp > 5: ro5_violations += 1
        if hba > 10: ro5_violations += 1
        if hbd > 5: ro5_violations += 1
        
        # Calculate QED (Quantitative Estimate of Drug-likeness)
        qed = QED.qed(mol)
        
        # Get formula
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        
        # Get Murcko Scaffold
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        
        # Calculate additional descriptors
        aromatic_proportion = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) / mol.GetNumHeavyAtoms() if mol.GetNumHeavyAtoms() > 0 else 0
        
        return {
            "smiles": smiles,
            "molecular_weight": mw,
            "logp": logp,
            "tpsa": tpsa,
            "num_h_acceptors": hba,
            "num_h_donors": hbd,
            "num_rotatable_bonds": rotatable_bonds,
            "ring_count": rings,
            "aromatic_rings": aromatic_rings,
            "heavy_atoms": heavy_atoms,
            "ro5_violations": ro5_violations,
            "qed": qed,
            "formula": formula,
            "scaffold": scaffold,
            "aromatic_proportion": aromatic_proportion,
            "has_3d_coords": has_3d
        }
    except Exception as e:
        return {"error": str(e)}

def main():
    """Process SMILES strings from command line arguments."""
    if len(sys.argv) < 2:
        print(json.dumps({"error": "No SMILES provided. Usage: python descriptors.py SMILES"}))
        sys.exit(1)
    
    smiles = sys.argv[1]
    result = calculate_descriptors(smiles)
    print(json.dumps(result))

if __name__ == "__main__":
    main() 