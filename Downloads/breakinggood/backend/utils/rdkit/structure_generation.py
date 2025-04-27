#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import base64
from io import BytesIO

def generate_3d_structure(smiles):
    """Generate 3D coordinates for a molecule from SMILES"""
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        success = AllChem.EmbedMolecule(mol, useRandomCoords=True)
        if success == -1:
            return {'error': 'Failed to generate 3D coordinates'}
        
        # Perform MMFF optimization
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Get atom positions
        conf = mol.GetConformer()
        positions = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            atom = mol.GetAtomWithIdx(i)
            element = atom.GetSymbol()
            positions.append({
                'atom_id': i,
                'element': element,
                'x': pos.x,
                'y': pos.y,
                'z': pos.z
            })
        
        # Calculate energy (MMFF94)
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
            energy = ff.CalcEnergy() if ff else None
        except:
            energy = None
        
        # Get bonds
        bonds = []
        for bond in mol.GetBonds():
            bond_type = str(bond.GetBondType())
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            bonds.append({
                'begin_atom': begin_atom,
                'end_atom': end_atom,
                'bond_type': bond_type
            })
        
        # Get 3D representation as PDB string
        pdb_string = Chem.MolToPDBBlock(mol)
        
        # Get 3D representation as XYZ format
        xyz_block = """{}\n\n""".format(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            atom = mol.GetAtomWithIdx(i)
            element = atom.GetSymbol()
            xyz_block += "{}\t{:.6f}\t{:.6f}\t{:.6f}\n".format(
                element, pos.x, pos.y, pos.z
            )
        
        # Get MOL block representation for 3Dmol.js visualization
        mol_block = Chem.MolToMolBlock(mol)
        
        return {
            'success': True,
            'atom_count': mol.GetNumAtoms(),
            'positions': positions,
            'bonds': bonds,
            'energy': energy,
            'pdb_block': pdb_string,
            'xyz_block': xyz_block,
            'molblock': mol_block  # Added for 3D visualization
        }
    
    except Exception as e:
        return {'error': str(e)}

def main():
    """Main function to process SMILES from command line arguments"""
    if len(sys.argv) < 2:
        print(json.dumps({"error": "No SMILES string provided"}))
        return
    
    smiles = sys.argv[1]
    result = generate_3d_structure(smiles)
    print(json.dumps(result))

if __name__ == '__main__':
    main()
