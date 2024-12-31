import numpy as np
import pandas as pd
import h5py

# File paths
METADATA_FILE = (
    "/projects/synsight/data/website_data/jump_compounds_matrix_metadata.parquet"
)
H5_DISTANCE_FILE = "/projects/synsight/data/website_data/nearest_neighbors.h5"

# Variables to hold the data
metadata = None
h5_file = None

def load_data():
    """Loads metadata and HDF5 file into memory if not already loaded."""
    global metadata, h5_file
    if metadata is None:
        metadata = pd.read_parquet(METADATA_FILE)
    if h5_file is None:
        h5_file = h5py.File(H5_DISTANCE_FILE, 'r')

def find_closest_molecules(query, n, query_by="id"):
    """Find the n closest molecules to the given query."""
    load_data()

    if query_by == "inchi":
        if query not in metadata["Metadata_InChI"].values:
            return {"error": "Query InChI not found"}
        query_id = metadata.loc[metadata["Metadata_InChI"] == query, "Metadata_JCP2022"].values[0]
    else:
        if query not in metadata["Metadata_JCP2022"].values:
            return {"error": "Query ID not found"}
        query_id = query

    if query_id not in h5_file:
        return {"error": f"Query ID {query_id} not found in distance data"}

    # Load the closest distances and molecule IDs from the HDF5 file
    closest_ids = h5_file[f"{query_id}/closest_ids"][:n+1].astype(str)
    closest_distances = h5_file[f"{query_id}/distances"][:n+1]

    # Retrieve metadata for the closest molecules
    closest_molecules = metadata[metadata["Metadata_JCP2022"].isin(closest_ids)].copy()
    closest_molecules.loc[:, "distance"] = [
        closest_distances[np.where(closest_ids == molecule_id)[0][0]]
        for molecule_id in closest_molecules["Metadata_JCP2022"].values
    ]

    return closest_molecules[
        ["Metadata_JCP2022", "Metadata_InChI", "distance"]
    ].sort_values("distance").to_dict(orient="records")

def find_distance_to_dmso(query, query_type):
    """Find the distance of the query molecule to DMSO."""
    load_data()

    if query_type == "inchi":
        if query not in metadata["Metadata_InChI"].values:
            return {"error": "Query InChI not found"}
        query_id = metadata.loc[metadata["Metadata_InChI"] == query, "Metadata_JCP2022"].values[0]
    else:
        if query not in metadata["Metadata_JCP2022"].values:
            return {"error": "Query ID not found"}
        query_id = query

    if query_id not in h5_file:
        return {"error": f"Query ID {query_id} not found in distance data"}

    # Load the distance to DMSO
    dmso_distance = h5_file[f"{query_id}/dmso_distance"][()]
    return float(dmso_distance)

# Ensure the HDF5 file is properly closed when the script finishes
import atexit
atexit.register(lambda: h5_file.close() if h5_file else None)
