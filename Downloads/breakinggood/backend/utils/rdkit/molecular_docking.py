#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import os
import tempfile
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def prepare_ligand(smiles):
    """Prepare a ligand from SMILES for docking"""
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, "Invalid SMILES string"
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        status = AllChem.EmbedMolecule(mol, randomSeed=42)
        if status == -1:
            return None, "Could not generate 3D coordinates"
        
        # Energy minimize the molecule
        AllChem.MMFFOptimizeMolecule(mol)
        
        return mol, None
    except Exception as e:
        return None, str(e)

def perform_docking(receptor_pdb, ligand_smiles, center=None, size=None, exhaustiveness=8):
    """Perform molecular docking using AutoDock Vina"""
    try:
        # Create temporary directory for docking files
        temp_dir = tempfile.mkdtemp()
        
        # Prepare ligand
        mol, error = prepare_ligand(ligand_smiles)
        if error:
            return {"error": error}
        
        # Write receptor to PDB file
        receptor_path = os.path.join(temp_dir, "receptor.pdb")
        with open(receptor_path, "w") as f:
            f.write(receptor_pdb)
        
        # Write ligand to PDBQT file
        ligand_pdb_path = os.path.join(temp_dir, "ligand.pdb")
        ligand_pdbqt_path = os.path.join(temp_dir, "ligand.pdbqt")
        Chem.MolToPDBFile(mol, ligand_pdb_path)
        
        # Convert PDB to PDBQT for receptor using MGLTools
        receptor_pdbqt_path = os.path.join(temp_dir, "receptor.pdbqt")
        os.system(f"python -m meeko.pdbqt receptor -r {receptor_path} -o {receptor_pdbqt_path}")
        
        # Convert PDB to PDBQT for ligand using MGLTools
        os.system(f"python -m meeko.pdbqt protein -i {ligand_pdb_path} -o {ligand_pdbqt_path}")
        
        # Determine center and size if not provided
        if center is None or size is None:
            # Use receptor center and dimensions
            # This is a simplified approach - in real applications, 
            # you would use a more sophisticated binding site detection
            coords = []
            with open(receptor_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM"):
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
            
            coords = np.array(coords)
            center = coords.mean(axis=0)
            size = coords.max(axis=0) - coords.min(axis=0) + 10  # Add 10 Angstrom padding
        
        # Create Vina configuration file
        config_path = os.path.join(temp_dir, "conf.txt")
        with open(config_path, "w") as f:
            f.write(f"receptor = {receptor_pdbqt_path}\n")
            f.write(f"ligand = {ligand_pdbqt_path}\n")
            f.write(f"center_x = {center[0]}\n")
            f.write(f"center_y = {center[1]}\n")
            f.write(f"center_z = {center[2]}\n")
            f.write(f"size_x = {size[0]}\n")
            f.write(f"size_y = {size[1]}\n")
            f.write(f"size_z = {size[2]}\n")
            f.write(f"exhaustiveness = {exhaustiveness}\n")
            f.write(f"num_modes = 9\n")
        
        # Output file for docking results
        output_path = os.path.join(temp_dir, "docking_out.pdbqt")
        
        # Run Vina docking
        cmd = f"vina --config {config_path} --out {output_path}"
        os.system(cmd)
        
        # Parse docking results
        results = []
        current_pose = {}
        with open(output_path, "r") as f:
            for line in f:
                if line.startswith("MODEL"):
                    current_pose = {"atoms": []}
                elif line.startswith("ENDMDL"):
                    results.append(current_pose)
                elif line.startswith("REMARK VINA RESULT"):
                    parts = line.split()
                    current_pose["affinity"] = float(parts[3])
                    current_pose["rmsd_lb"] = float(parts[4])
                    current_pose["rmsd_ub"] = float(parts[5])
                elif line.startswith("ATOM") or line.startswith("HETATM"):
                    atom = {
                        "serial": int(line[6:11]),
                        "name": line[12:16].strip(),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "element": line[76:78].strip()
                    }
                    current_pose["atoms"].append(atom)
        
        # Create a simplified response for best pose
        if results:
            best_pose = results[0]  # First pose is the best one
            response = {
                "smiles": ligand_smiles,
                "binding_affinity": best_pose["affinity"],
                "rmsd_lb": best_pose["rmsd_lb"],
                "rmsd_ub": best_pose["rmsd_ub"],
                "atom_count": len(best_pose["atoms"]),
                "poses": [{
                    "affinity": pose["affinity"],
                    "rmsd_lb": pose["rmsd_lb"],
                    "rmsd_ub": pose["rmsd_ub"]
                } for pose in results],
                "success": True
            }
        else:
            response = {"error": "No docking results obtained", "success": False}
        
        # Clean up temp directory
        shutil.rmtree(temp_dir)
        
        return response
    except Exception as e:
        # Clean up temp directory if it exists
        if 'temp_dir' in locals():
            shutil.rmtree(temp_dir)
        return {"error": str(e), "success": False}

def main():
    """Main function to process docking inputs from command line arguments"""
    if len(sys.argv) < 3:
        print(json.dumps({"error": "Need receptor file path and ligand SMILES"}))
        sys.exit(1)
    
    receptor_path = sys.argv[1]
    ligand_smiles = sys.argv[2]
    
    # Optional parameters
    center = None
    size = None
    exhaustiveness = 8
    
    if len(sys.argv) > 5:
        center = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
    
    if len(sys.argv) > 8:
        size = [float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])]
    
    if len(sys.argv) > 9:
        exhaustiveness = int(sys.argv[9])
    
    # Read receptor file
    with open(receptor_path, "r") as f:
        receptor_pdb = f.read()
    
    # Perform docking
    results = perform_docking(receptor_pdb, ligand_smiles, center, size, exhaustiveness)
    print(json.dumps(results))

if __name__ == "__main__":
    main()
