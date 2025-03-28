#!/usr/bin/env python
import pandas as pd
from rdkit import Chem
import pubchempy as pcp
from joblib import Parallel, delayed
from tqdm import tqdm
import time
import random

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
    """
    Lookup PubChem compound ID using PubChemPy.
    If the PubChem server returns a 'PUGREST.ServerBusy' error, wait with an added random jitter before retrying.
    """
    max_retries = 7
    delay = 0.5  # initial delay in seconds
    for attempt in range(max_retries):
        try:
            compounds = pcp.get_compounds(inchi_str, 'inchi')
            if compounds:
                return compounds[0].cid
            return None
        except Exception as e:
            if "PUGREST.ServerBusy" in str(e):
                # Add a random jitter to avoid simultaneous retries across workers
                jitter = random.uniform(0, 1)
                sleep_time = delay + 10*jitter
             #   print(f"PubChem server busy, attempt {attempt+1}/{max_retries}. Retrying in {sleep_time:.2f} seconds...")
                time.sleep(sleep_time)
                delay *= 2  # exponential backoff
            else:
                print(f"Error fetching PubChem ID: {e}")
                return None
    return None

def get_chembl_id(inchikey):
    """
    Lookup ChEMBL ID using the ChEMBL REST API.
    Import new_client inside the function to avoid pickling issues.
    """
    try:
        from chembl_webresource_client.new_client import new_client
        results = new_client.molecule.filter(molecule_structures__standard_inchi_key=inchikey)
        if results:
            return results[0]['molecule_chembl_id']
    except Exception as e:
        print(f"Error fetching ChEMBL ID: {e}")
    return None

def process_molecule(inchi_str):
    """Process one molecule: convert InChI and lookup PubChem and ChEMBL IDs."""
    ikey = get_inchikey(inchi_str)
    if ikey:
        pub_id = get_pubchem_id(inchi_str)
        chembl_id = get_chembl_id(ikey)
    else:
        pub_id, chembl_id = None, None
    return ikey, pub_id, chembl_id

def add_pubchem_chembl_ids(df, use_threading=True):
    """
    Add PubChem and ChEMBL IDs to the DataFrame by processing each molecule in parallel.
    
    use_threading: if True, uses the threading backend to avoid pickling issues.
    """
    backend = "threading" if use_threading else "loky"
    results = Parallel(n_jobs=10, backend=backend)(
        delayed(process_molecule)(inchi_str) for inchi_str in tqdm(df['Metadata_InChI'], desc="Processing molecules")
    )
    inchikeys, pubchem_ids, chembl_ids = zip(*results)
    df['InChIKey'] = inchikeys
    df['PubChem_ID'] = pubchem_ids
    df['ChEMBL_ID'] = chembl_ids
    return df

if __name__ == "__main__":
    METADATA_FILE = ("/projects/synsight/data/website_data/jump_compounds_matrix_metadata.parquet")
    df = pd.read_parquet(METADATA_FILE)
    
    updated_df = add_pubchem_chembl_ids(df, use_threading=True)
    updated_df.to_parquet("compounds_with_pubchem_chembl.parquet", index=False)
    print("Updated DataFrame saved to molecules_with_pubchem_chembl.csv")
