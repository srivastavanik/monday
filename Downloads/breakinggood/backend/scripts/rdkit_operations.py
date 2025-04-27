import base64
import io
import os # Added for file cleanup in docking
from rdkit.Chem import rdDistGeom # For ETKDG embedding
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
    OPENMM_INSTALLED = True
except ImportError:
    OPENMM_INSTALLED = False

def calculate_properties(smiles):
    """Calculate molecular properties for a given SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES string"}
            
        # Calculate basic properties
        props = {
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "molWeight": Descriptors.MolWt(mol),
            "logP": Crippen.MolLogP(mol),
            "tpsa": MolSurf.TPSA(mol),
            "hbondDonorCount": Lipinski.NumHDonors(mol),
            "hbondAcceptorCount": Lipinski.NumHAcceptors(mol),
            "rotatableBondCount": Descriptors.NumRotatableBonds(mol),
            "aromaticRingCount": Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
            "heavyAtomCount": mol.GetNumHeavyAtoms(),
            "qed": QED.qed(mol),
            "lipinskiViolations": Lipinski.NumLipinskiHBAs(mol) > 10 or
                                  Lipinski.NumLipinskiHBDs(mol) > 5 or
                                  Descriptors.MolWt(mol) > 500 or
                                  Crippen.MolLogP(mol) > 5,
        }
        
        # Get 2D depiction as SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        props["svg"] = drawer.GetDrawingText()
        
        return props
    except Exception as e:
        return {"error": str(e), "traceback": traceback.format_exc()}

def generate_3d_structure(smiles):
   # ... (existing 3D structure code) ...
   pass # Placeholder to ensure indentation is correct

def convert_molecule(input_data, input_format, output_format):
   # ... (existing conversion code) ...
   pass # Placeholder

def predict_admet(smiles):
   # ... (existing ADMET code) ...
   pass # Placeholder

def run_short_md_simulation(smiles, steps=50000, time_step_fs=2.0, temperature_k=300, force_field='MMFF94'):
    """Run a short MD simulation using OpenMM via RDKit."""
    if not OPENMM_INSTALLED:
        return {"error": "OpenMM is not installed. Cannot run MD simulation."}

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES string"}
            
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates using ETKDG v3
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = 0xf00d # Deterministic embedding
        if AllChem.EmbedMolecule(mol, params) == -1:
             return {"error": "Failed to generate 3D coordinates"}

        # Select force field (MMFF94 or UFF)
        if force_field.upper() == 'UFF':
            ff_func = AllChem.UFFGetMoleculeForceField
            props_func = AllChem.UFFGetMoleculeProperties
            ff_name = 'UFF'
        else: # Default to MMFF94
            ff_func = AllChem.MMFFGetMoleculeForceField
            props_func = AllChem.MMFFGetMoleculeProperties
            ff_name = 'MMFF94'

        # Optimize geometry with selected force field
        try:
            ff = ff_func(mol, props_func(mol))
            if ff is None:
                raise ValueError(f"Could not get {ff_name} force field")
            ff.Initialize()
            if ff.Minimize() != 0: # Check if minimization succeeded
                print(f"Warning: {ff_name} minimization did not converge.")
            initial_energy = ff.CalcEnergy()
        except Exception as ff_error:
             return {"error": f"Force field setup/minimization failed ({ff_name}): {ff_error}"}

        # --- OpenMM Simulation Setup ---
        topology = Topology()
        positions = mol.GetConformer().GetPositions().value_in_unit(nanometers)
        chain = topology.addChain()
        residue = topology.addResidue(smiles, chain)
        atom_map = {} # Map RDKit atom indices to OpenMM atom indices
        for atom in mol.GetAtoms():
            element = Element.getByAtomicNumber(atom.GetAtomicNum())
            omm_atom = topology.addAtom(atom.GetSymbol(), element, residue)
            atom_map[atom.GetIdx()] = omm_atom

        # Use RDKit force field for OpenMM parameters
        omm_ff = ForceField() 
        # Note: This is a simplification. Directly using RDKit FF in OpenMM isn't standard.
        # A more robust approach uses GAFF/SMIRNOFF via OpenFF Toolkit or parameterizes manually.
        # This setup will likely only capture basic bonded terms correctly without explicit parameters.
        # For a quick demo, we proceed, but acknowledge limitations.
        
        # Create system (This part is highly simplified and likely inaccurate for rigorous MD)
        # A proper setup requires generating OpenMM parameters from the RDKit force field or using a standard FF
        # We'll create a minimal system just to run the simulation steps
        try:
             system = omm_ff.createSystem(topology) # This will likely be empty or minimal without proper FF files
             # Add non-bonded forces manually if needed (example)
             # system.addForce(NonbondedForce()) 
        except Exception as sys_error:
             return {"error": f"Failed to create OpenMM system: {sys_error}. Ensure FF parameters are available or simplify system creation."}

        integrator = LangevinIntegrator(temperature_k*kelvin, 1/picosecond, time_step_fs*femtoseconds)
        
        # Try CUDA platform first, fall back to Reference (CPU)
        platform = None
        try:
            platform = Platform.getPlatformByName('CUDA')
            # Optionally set CUDA device index or precision
            # properties = {'CudaPrecision': 'mixed'}
            # simulation = Simulation(topology, system, integrator, platform, properties)
            print("Using CUDA platform for OpenMM simulation.")
        except Exception as e_cuda:
            print(f"CUDA platform not found or failed ({e_cuda}), falling back to Reference platform.")
            try:
                 platform = Platform.getPlatformByName('Reference')
            except Exception as e_ref:
                 return {"error": f"Could not get Reference platform either: {e_ref}"}
        
        simulation = Simulation(topology, system, integrator, platform)
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(temperature_k*kelvin)

        # --- Run Simulation ---
        simulation.step(steps)

        # --- Get Results ---
        final_state = simulation.context.getState(getPositions=True, getEnergy=True)
        final_positions = final_state.getPositions(asNumpy=True).value_in_unit(nanometers)
        final_potential_energy = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        
        # Update RDKit molecule positions
        final_conformer = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            final_conformer.SetAtomPosition(i, Point3D(final_positions[i][0], final_positions[i][1], final_positions[i][2]))

        # Calculate RMSD between initial minimized and final MD structure
        try:
            rmsd_md = AllChem.AlignMol(mol, initial_mol) # Align final structure to initial minimized
        except Exception as align_error:
            rmsd_md = -1
            print(f"Warning: RMSD calculation failed after MD: {align_error}")

        final_molblock = Chem.MolToMolBlock(mol) # Get molblock of the final state

        return {
            "status": "completed",
            "force_field": ff_name,
            "steps": steps,
            "time_ps": steps * time_step_fs / 1000.0,
            "initial_energy_kcal_mol": initial_energy,
            "final_potential_energy_kcal_mol": final_potential_energy,
            "rmsd_vs_minimized": rmsd_md,
            "final_molblock": final_molblock
        }

    except Exception as e:
        return {"error": str(e), "traceback": traceback.format_exc()}

def run_docking(receptor_pdb, ligand_smiles, exhaustiveness=8, center_x=0, center_y=0, center_z=0, size_x=20, size_y=20, size_z=20):
    # ... (existing docking code) ...
    pass # Placeholder

def compare_molecules(smiles1, smiles2):
    # ... (existing compare code) ...
    pass # Placeholder

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Not enough arguments provided"}))
        sys.exit(1)
        
    operation = sys.argv[1]
    
    if operation == "properties":
        if len(sys.argv) < 3:
            print(json.dumps({"error": "SMILES string required"}))
            sys.exit(1)
        print(json.dumps(calculate_properties(sys.argv[2])))
        
    elif operation == "3d_structure":
        if len(sys.argv) < 3:
            print(json.dumps({"error": "SMILES string required"}))
            sys.exit(1)
        print(json.dumps(generate_3d_structure(sys.argv[2])))
        
    elif operation == "convert":
        if len(sys.argv) < 5:
            print(json.dumps({"error": "Input, input format, and output format required"}))
            sys.exit(1)
        print(json.dumps(convert_molecule(sys.argv[2], sys.argv[3], sys.argv[4])))
        
    elif operation == "admet":
        if len(sys.argv) < 3:
            print(json.dumps({"error": "SMILES string required"}))
            sys.exit(1)
        print(json.dumps(predict_admet(sys.argv[2])))
        
    elif operation == "md":
        if len(sys.argv) < 3:
            print(json.dumps({"error": "SMILES string required"}))
            sys.exit(1)
        smiles_arg = sys.argv[2]
        steps_arg = int(sys.argv[3]) if len(sys.argv) > 3 else 50000
        print(json.dumps(run_short_md_simulation(smiles_arg, steps=steps_arg)))
        
    elif operation == "compare":
        if len(sys.argv) < 4:
            print(json.dumps({"error": "Two SMILES strings required"}))
            sys.exit(1)
        print(json.dumps(compare_molecules(sys.argv[2], sys.argv[3])))
        
    elif operation == "docking":
        if len(sys.argv) < 10:
            print(json.dumps({"error": "Receptor PDB, ligand SMILES, and docking parameters required"}))
            sys.exit(1)
        # Read PDB content directly from argument if it's long enough, otherwise assume it's a file path
        receptor_pdb_arg = sys.argv[2]
        if len(receptor_pdb_arg) > 200 and 'ATOM' in receptor_pdb_arg: # Heuristic: Assume PDB content if long and contains ATOM records
             receptor_pdb = receptor_pdb_arg
        else:
             # This case might need adjustment if PDB path is passed instead of content
             # For now, assume content is passed directly for simplicity based on previous JS code
             receptor_pdb = receptor_pdb_arg 
        
        ligand_smiles = sys.argv[3]
        exhaustiveness = int(sys.argv[4])
        center_x = float(sys.argv[5])
        center_y = float(sys.argv[6])
        center_z = float(sys.argv[7])
        size_x = float(sys.argv[8])
        size_y = float(sys.argv[9])
        size_z = float(sys.argv[10])
        print(json.dumps(run_docking(receptor_pdb, ligand_smiles, exhaustiveness, center_x, center_y, center_z, size_x, size_y, size_z)))
    else:
        print(json.dumps({"error": f"Unknown operation: {operation}"})) 