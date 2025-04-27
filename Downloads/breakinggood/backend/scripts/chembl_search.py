
import sys
import json
import traceback
from chembl_webresource_client.new_client import new_client

def search_chembl(query_type, query_value, limit=10):
    try:
        results = []
        
        if query_type == 'molecule':
            molecule = new_client.molecule
            if query_value.startswith('CHEMBL'):
                res = molecule.filter(chembl_id=query_value).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
            else:
                res = molecule.filter(pref_name__icontains=query_value).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
            results = list(res)[:limit]
            
        elif query_type == 'target':
            target = new_client.target
            if query_value.startswith('CHEMBL'):
                res = target.filter(target_chembl_id=query_value)
            else:
                res = target.filter(pref_name__icontains=query_value)
            results = list(res)[:limit]
            
        elif query_type == 'activity':
            activity = new_client.activity
            if query_value.startswith('CHEMBL'):
                res = activity.filter(molecule_chembl_id=query_value).only(['molecule_chembl_id', 'target_chembl_id', 'standard_type', 'standard_value', 'standard_units'])
            else:
                res = activity.filter(pref_name__icontains=query_value).only(['molecule_chembl_id', 'target_chembl_id', 'standard_type', 'standard_value', 'standard_units'])
            results = list(res)[:limit]
            
        elif query_type == 'similarity':
            similarity = new_client.similarity
            res = similarity.filter(smiles=query_value, similarity=85).only(['molecule_chembl_id', 'pref_name', 'molecule_structures', 'similarity'])
            results = list(res)[:limit]
            
        elif query_type == 'substructure':
            substructure = new_client.substructure
            res = substructure.filter(smiles=query_value).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])
            results = list(res)[:limit]
            
        return json.dumps(results)
    except Exception as e:
        return json.dumps({"error": str(e), "traceback": traceback.format_exc()})

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(json.dumps({"error": "Not enough arguments provided"}))
        sys.exit(1)
        
    query_type = sys.argv[1]
    query_value = sys.argv[2]
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    
    print(search_chembl(query_type, query_value, limit))
