#!/usr/bin/env python3
import sys
import json
import argparse
from rdkit import Chem

def convert_format(input_data, input_format, output_format):
    """Convert between molecular formats (SMILES, Molfile, InChI)."""
    mol = None
    
    # Read molecule based on input format
    if input_format == 'smiles':
        mol = Chem.MolFromSmiles(input_data)
    elif input_format == 'mol':
        mol = Chem.MolFromMolBlock(input_data, removeHs=False) # Keep Hs for Molfile input
    elif input_format == 'inchi':
        mol = Chem.MolFromInchi(input_data)
    else:
        return {"error": f"Unsupported input format: {input_format}"}

    if mol is None:
        return {"error": f"Failed to parse input molecule with format '{input_format}'"}

    # Generate output based on requested format
    output_data = None
    if output_format == 'smiles':
        # Generate canonical SMILES
        output_data = Chem.MolToSmiles(mol, isomericSmiles=True)
    elif output_format == 'mol':
        # Add explicit hydrogens if needed before generating Molfile
        # mol = Chem.AddHs(mol)
        output_data = Chem.MolToMolBlock(mol)
    elif output_format == 'inchi':
        output_data = Chem.MolToInchi(mol)
    elif output_format == 'inchikey':
        output_data = Chem.MolToInchiKey(mol)
    else:
        return {"error": f"Unsupported output format: {output_format}"}

    if output_data is None:
         return {"error": f"Failed to generate output format '{output_format}'"}
         
    return {"output": output_data}

def main():
    parser = argparse.ArgumentParser(description='RDKit Molecular Format Converter')
    parser.add_argument('input_data', type=str, help='Input molecule data (SMILES, Molfile block, or InChI)')
    parser.add_argument('input_format', type=str, choices=['smiles', 'mol', 'inchi'],
                        help='Format of the input data')
    parser.add_argument('output_format', type=str, choices=['smiles', 'mol', 'inchi', 'inchikey'],
                        help='Desired output format')
    
    args = parser.parse_args()
    
    result = convert_format(args.input_data, args.input_format, args.output_format)
    print(json.dumps(result))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1) 