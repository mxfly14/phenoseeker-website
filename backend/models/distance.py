import numpy as np
import pandas as pd

# Load data
METADATA_FILE = "/projects/synsight/data/website_data/jump_compounds_matrix_metadata.parquet"
DISTANCE_MATRIX_FILE = "/projects/synsight/data/website_data/jump_compounds_matrix.npy"

metadata = pd.read_parquet(METADATA_FILE)
distance_matrix = np.load(DISTANCE_MATRIX_FILE)

def find_closest_molecules(pos_control_id, n):
    # Find the index of the positive control
    if pos_control_id not in metadata['Metadata_JCP2022'].values:
        return {'error': 'Positive control ID not found'}
    
    pos_control_idx = metadata.index[metadata['id'] == pos_control_id][0]
    
    # Get distances for the positive control
    distances = distance_matrix[pos_control_idx]
    
    # Find the indices of the n smallest distances (excluding the control itself)
    closest_indices = np.argsort(distances)[1:n + 1]
    
    # Retrieve molecule IDs and distances
    closest_molecules = metadata.iloc[closest_indices]
    closest_molecules['distance'] = distances[closest_indices]
    
    return closest_molecules[['Metadata_JCP2022', 'Metadata_InChI', 'distance']].to_dict(orient='records')
