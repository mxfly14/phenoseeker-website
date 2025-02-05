import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

DATA_NPY_FILE = "/projects/synsight/data/website_data/jump_compounds_matrix.npy"
METADATA_FILE = "/projects/synsight/data/website_data/jump_compounds_matrix_metadata.parquet"
output_file = '/projects/synsight/data/website_data/nearest_neighbors.h5'  # Output HDF5 file
m = 1000  # Number of nearest neighbors


# Charger la matrice depuis le fichier .npy
matrix = np.load(DATA_NPY_FILE)
all_values = matrix.flatten()

# Load metadata
metadata = pd.read_parquet(METADATA_FILE)
metadata_ids = metadata['Metadata_JCP2022'].values  # Unique molecule IDs

def process_row(i):
    """
    Process a single row of the distance matrix to find the m closest neighbors
    and include the distance to a specific molecule.

    Args:
        i (int): Index of the row in the distance matrix.

    Returns:
        tuple: (molecule_id, closest_ids, closest_distances)
    """
    distances = matrix[i]

    # Distance to the target molecule (JCP2022_033924)
    try:
        target_index = np.where(metadata_ids == 'JCP2022_033924')[0][0]
    except IndexError:
        raise ValueError("Molecule JCP2022_033924 (DMSO) not found in metadata")
    dmso_distance = distances[target_index]

    # Find m closest molecules using partial sorting (excluding self if needed)
    closest_indices = np.argpartition(distances, m)[1:m+1]  # Top m indices (unsorted)
    closest_distances = distances[closest_indices]

    # Sort these m indices to ensure proper order
    sorted_indices_within_chunk = np.argsort(closest_distances)
    closest_indices = closest_indices[sorted_indices_within_chunk]
    closest_distances = closest_distances[sorted_indices_within_chunk]

    # Get IDs for the closest molecules
    closest_ids = metadata_ids[closest_indices]

    # Return the results
    return metadata_ids[i], closest_ids, closest_distances, dmso_distance

# Parallel processing
with Pool(processes=cpu_count()) as pool:
    # Use tqdm for progress tracking
    results = list(tqdm(pool.imap(process_row, range(matrix.shape[0])), total=matrix.shape[0]))


# Save results to HDF5
with h5py.File(output_file, 'w') as h5f:
    for molecule_id, closest_ids, closest_distances, dmso_distance in results:
        group = h5f.create_group(molecule_id)
        group.create_dataset('closest_ids', data=closest_ids.astype('S'))  # Save IDs as strings
        group.create_dataset('distances', data=closest_distances)
        group.create_dataset('dmso_distance', data=dmso_distance)