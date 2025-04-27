#!/usr/bin/env python3
import json
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

class MolecularSimilaritySearch:
    """A class for performing molecular similarity searches with various algorithms and fingerprints"""
    
    FINGERPRINT_TYPES = {
        'morgan': 'Morgan (ECFP4) fingerprint',
        'maccs': 'MACCS keys',
        'rdkit': 'RDKit topological fingerprint',
        'avalon': 'Avalon fingerprint'
    }
    
    SIMILARITY_METRICS = {
        'tanimoto': 'Tanimoto coefficient (Jaccard index)',
        'dice': 'Dice coefficient',
        'cosine': 'Cosine similarity',
        'sokal': 'Sokal similarity',
        'russel': 'Russel similarity'
    }
    
    def __init__(self, fingerprint_type='morgan', similarity_metric='tanimoto'):
        """Initialize the similarity search with specific fingerprint type and similarity metric"""
        self.set_fingerprint_type(fingerprint_type)
        self.set_similarity_metric(similarity_metric)
    
    def set_fingerprint_type(self, fingerprint_type):
        """Set the fingerprint type to use for similarity calculations"""
        if fingerprint_type not in self.FINGERPRINT_TYPES:
            raise ValueError(f"Unsupported fingerprint type: {fingerprint_type}. " 
                            f"Choose from: {', '.join(self.FINGERPRINT_TYPES.keys())}")
        self.fingerprint_type = fingerprint_type
    
    def set_similarity_metric(self, similarity_metric):
        """Set the similarity metric to use for comparisons"""
        if similarity_metric not in self.SIMILARITY_METRICS:
            raise ValueError(f"Unsupported similarity metric: {similarity_metric}. "
                            f"Choose from: {', '.join(self.SIMILARITY_METRICS.keys())}")
        self.similarity_metric = similarity_metric
    
    def generate_fingerprint(self, mol):
        """Generate a fingerprint based on the selected fingerprint type"""
        if mol is None:
            return None
            
        if self.fingerprint_type == 'morgan':
            return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        elif self.fingerprint_type == 'maccs':
            return MACCSkeys.GenMACCSKeys(mol)
        elif self.fingerprint_type == 'rdkit':
            return FingerprintMols.FingerprintMol(mol)
        elif self.fingerprint_type == 'avalon':
            try:
                from rdkit.Avalon import pyAvalonTools
                return pyAvalonTools.GetAvalonFP(mol)
            except ImportError:
                raise ImportError("Avalon fingerprints require the rdkit.Avalon module")
    
    def calculate_similarity(self, fp1, fp2):
        """Calculate similarity between two fingerprints using the selected metric"""
        if fp1 is None or fp2 is None:
            return 0.0
            
        if self.similarity_metric == 'tanimoto':
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        elif self.similarity_metric == 'dice':
            return DataStructs.DiceSimilarity(fp1, fp2)
        elif self.similarity_metric == 'cosine':
            return DataStructs.CosineSimilarity(fp1, fp2)
        elif self.similarity_metric == 'sokal':
            return DataStructs.SokalSimilarity(fp1, fp2)
        elif self.similarity_metric == 'russel':
            return DataStructs.RusselSimilarity(fp1, fp2)
    
    def search_by_similarity(self, query_smiles, target_library, threshold=0.7, max_results=50):
        """
        Search a library of compounds for similarity to a query molecule
        
        Args:
            query_smiles (str): SMILES string of query molecule
            target_library (list): List of dictionaries with 'smiles' and other properties
            threshold (float): Minimum similarity threshold (0.0-1.0)
            max_results (int): Maximum number of results to return
            
        Returns:
            list: Sorted list of similar compounds with similarity scores
        """
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            return {"error": "Invalid query SMILES string"}
            
        query_fp = self.generate_fingerprint(query_mol)
        results = []
        
        for compound in target_library:
            if 'smiles' not in compound:
                continue
                
            target_mol = Chem.MolFromSmiles(compound['smiles'])
            if target_mol is None:
                continue
                
            target_fp = self.generate_fingerprint(target_mol)
            similarity = self.calculate_similarity(query_fp, target_fp)
            
            if similarity >= threshold:
                result = compound.copy()
                result['similarity'] = similarity
                results.append(result)
        
        # Sort by similarity score (descending)
        results.sort(key=lambda x: x['similarity'], reverse=True)
        
        return results[:max_results]
    
    def search_multiple_references(self, reference_smiles_list, target_library, 
                                  threshold=0.7, max_results=50, aggregation='max'):
        """
        Search using multiple reference compounds
        
        Args:
            reference_smiles_list (list): List of SMILES strings for reference compounds
            target_library (list): List of dictionaries with 'smiles' and other properties
            threshold (float): Minimum similarity threshold (0.0-1.0)
            max_results (int): Maximum number of results to return
            aggregation (str): How to aggregate scores - 'max', 'mean', or 'min'
            
        Returns:
            list: Sorted list of similar compounds with similarity scores
        """
        if not reference_smiles_list:
            return {"error": "No reference SMILES strings provided"}
            
        # Generate fingerprints for reference compounds
        reference_fps = []
        for smiles in reference_smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            reference_fps.append(self.generate_fingerprint(mol))
            
        if not reference_fps:
            return {"error": "No valid reference molecules"}
            
        results = []
        for compound in target_library:
            if 'smiles' not in compound:
                continue
                
            target_mol = Chem.MolFromSmiles(compound['smiles'])
            if target_mol is None:
                continue
                
            target_fp = self.generate_fingerprint(target_mol)
            
            # Calculate similarity to each reference
            similarities = [self.calculate_similarity(ref_fp, target_fp) for ref_fp in reference_fps]
            
            # Aggregate similarities
            if aggregation == 'max':
                final_similarity = max(similarities)
            elif aggregation == 'mean':
                final_similarity = sum(similarities) / len(similarities)
            elif aggregation == 'min':
                final_similarity = min(similarities)
            else:
                final_similarity = max(similarities)  # Default to max
                
            if final_similarity >= threshold:
                result = compound.copy()
                result['similarity'] = final_similarity
                result['individual_similarities'] = similarities
                results.append(result)
        
        # Sort by aggregated similarity score (descending)
        results.sort(key=lambda x: x['similarity'], reverse=True)
        
        return results[:max_results]
    
    def diversity_selection(self, smiles_list, num_picks=10):
        """
        Select a diverse subset of molecules using MaxMin algorithm
        
        Args:
            smiles_list (list): List of SMILES strings
            num_picks (int): Number of diverse compounds to select
            
        Returns:
            list: Indices of selected diverse compounds
        """
        if not smiles_list:
            return {"error": "Empty SMILES list provided"}
            
        mols = []
        valid_indices = []
        
        # Create valid molecules and track their indices
        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mols.append(mol)
                valid_indices.append(i)
                
        if not mols:
            return {"error": "No valid molecules in the provided list"}
            
        if len(mols) <= num_picks:
            return valid_indices
            
        # Generate fingerprints
        fps = [self.generate_fingerprint(mol) for mol in mols]
        
        # Distance function for MaxMinPicker
        def distance_func(i, j):
            return 1.0 - self.calculate_similarity(fps[i], fps[j])
            
        # Use MaxMinPicker for diversity selection
        picker = MaxMinPicker()
        pick_indices = picker.LazyPick(distance_func, len(fps), num_picks, seed=42)
        
        # Map picked indices back to original indices
        return [valid_indices[idx] for idx in pick_indices]


def parse_args():
    parser = argparse.ArgumentParser(description='RDKit Molecular Similarity Search')
    parser.add_argument('--query', type=str, help='SMILES string of query molecule')
    parser.add_argument('--target_file', type=str, help='JSON file with target compounds')
    parser.add_argument('--fingerprint', type=str, default='morgan', 
                        choices=['morgan', 'maccs', 'rdkit', 'avalon'],
                        help='Fingerprint type to use')
    parser.add_argument('--metric', type=str, default='tanimoto',
                        choices=['tanimoto', 'dice', 'cosine', 'sokal', 'russel'],
                        help='Similarity metric to use')
    parser.add_argument('--threshold', type=float, default=0.7,
                        help='Similarity threshold (0.0-1.0)')
    parser.add_argument('--max_results', type=int, default=50,
                        help='Maximum number of results to return')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load target compounds
    with open(args.target_file, 'r') as f:
        target_library = json.load(f)
    
    # Initialize similarity search
    sim_search = MolecularSimilaritySearch(
        fingerprint_type=args.fingerprint,
        similarity_metric=args.metric
    )
    
    # Perform search
    results = sim_search.search_by_similarity(
        args.query,
        target_library,
        threshold=args.threshold,
        max_results=args.max_results
    )
    
    # Output results
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main() 