#!/usr/bin/env python
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi
from tqdm import tqdm
import pubchempy as pcp
from chembl_webresource_client.new_client import new_client

def get_inchikey(inchi_str):
    """Convert InChI to InChIKey using RDKit."""
    try:
        mol = Chem.MolFromInchi(inchi_str)
        if mol is None:
            return None
        return Chem.MolToInchiKey(mol)
    except Exception as e:
        print(f"Error converting InChI to InChIKey: {e}")
        return None

def get_pubchem_id(inchi_str):
    """Lookup PubChem compound ID using PubChemPy."""
    try:
        compounds = pcp.get_compounds(inchi_str, 'inchi')
        if compounds:
            return compounds[0].cid
    except Exception as e:
        print(f"Error fetching PubChem ID: {e}")
    return None

def get_chembl_id(inchikey):
    """Lookup ChEMBL ID using the ChEMBL REST API."""
    try:
        results = new_client.molecule.filter(molecule_structures__standard_inchi_key=inchikey)
        if results:
            return results[0]['molecule_chembl_id']
    except Exception as e:
        print(f"Error fetching ChEMBL ID: {e}")
    return None

def add_pubchem_chembl_ids(df):
    """Add PubChem and ChEMBL IDs to the DataFrame."""
    # Lists to hold the results
    inchikeys = []
    pubchem_ids = []
    chembl_ids = []
    
    for idx, inchi_str in tqdm(enumerate(df['Metadata_InChI'])):
        ikey = get_inchikey(inchi_str)
        inchikeys.append(ikey)
        
        if ikey:
            pub_id = get_pubchem_id(inchi_str)
            chembl_id = get_chembl_id(ikey)
        else:
            pub_id, chembl_id = None, None
            
        pubchem_ids.append(pub_id)
        chembl_ids.append(chembl_id)
        
        if idx % 1000 == 0:
            print(f"Processed {idx} molecules")
    
    # Add the new columns to the DataFrame
    df['InChIKey'] = inchikeys
    df['PubChem_ID'] = pubchem_ids
    df['ChEMBL_ID'] = chembl_ids
    return df

if __name__ == "__main__":
    METADATA_FILE = ("/projects/synsight/data/website_data/jump_compounds_matrix_metadata.parquet")
    df = pd.read_parquet(METADATA_FILE)
        
    updated_df = add_pubchem_chembl_ids(df)
    updated_df.to_parquet("compounds_with_pubchem_chembl_slow.parquet", index=False)
    print("Updated DataFrame saved to molecules_with_pubchem_chembl.csv")
