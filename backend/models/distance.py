import numpy as np
import pandas as pd

# Load data
METADATA_FILE = (
    "/projects/synsight/data/website_data/jump_compounds_matrix_metadata.parquet"
)
DISTANCE_MATRIX_FILE = "/projects/synsight/data/website_data/jump_compounds_matrix.npy"
# Variables to hold the data
metadata = None
distance_matrix = None


def load_data():
    """Loads metadata and distance matrix into memory if not already loaded."""
    global metadata, distance_matrix
    if metadata is None:
        metadata = pd.read_parquet(METADATA_FILE)
    if distance_matrix is None:
        distance_matrix = np.load(DISTANCE_MATRIX_FILE)


def find_closest_molecules(query, n, query_by="id"):
    """Find the n closest molecules to the given query."""
    load_data()

    if query_by == "inchi":
        if query not in metadata["Metadata_InChI"].values:
            return {"error": "Query InChI not found"}
        query_idx = metadata.index[metadata["Metadata_InChI"] == query][0]
    else:
        if query not in metadata["Metadata_JCP2022"].values:
            return {"error": "Query ID not found"}
        query_idx = metadata.index[metadata["Metadata_JCP2022"] == query][0]

    # Get distances for the query molecule
    distances = distance_matrix[query_idx]

    # Find the indices of the n smallest distances (excluding the query itself)
    closest_indices = np.argsort(distances)[1 : n + 1]

    # Retrieve molecule IDs, InChIs, and distances
    closest_molecules = metadata.iloc[closest_indices]
    closest_molecules["distance"] = distances[closest_indices]

    return closest_molecules[
        ["Metadata_JCP2022", "Metadata_InChI", "distance"]
    ].to_dict(orient="records")


def find_distance_to_dmso(query, query_type):
    """Find the distance of the query molecule to DMSO."""
    load_data()

    if query_type == "inchi":
        query_idx = metadata.index[metadata["Metadata_InChI"] == query][0]
    else:
        query_idx = metadata.index[metadata["Metadata_JCP2022"] == query][0]

    # Find the index for DMSO (assuming "DMSO" is the ID for DMSO in the metadata)
    dmso_idx = metadata.index[metadata["Metadata_InChI"] == "JCP2022_033924"][0]

    return float(distance_matrix[query_idx, dmso_idx])
