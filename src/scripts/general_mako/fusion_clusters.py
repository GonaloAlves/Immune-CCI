# Import packages
import os
import scanpy as sc
import time

# This script will merge cells from clusters from the user agument dataset

# Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    print("Loading h5ad file")
    return sc.read_h5ad(file_path)


def merge_clusters(adata, resolution, merge_groups, new_resolution):
    """
    Merge multiple clusters into specified target clusters and save the result into a new resolution column.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    resolution (str): The original resolution column in `adata.obs` to use as a base.
    merge_groups (list of tuples): List of tuples where each tuple specifies:
                                   (list of clusters to merge, target cluster).
    new_resolution (str): The name of the new resolution column to store the result.

    Returns:
    AnnData: The modified AnnData object with the new resolution column.
    """
    print(f"Merging clusters in resolution '{resolution}'...")

    # Create a copy of the original resolution column
    adata.obs[new_resolution] = adata.obs[resolution]

    # Iterate through each group of clusters to merge
    for clusters_from, cluster_to in merge_groups:
        print(f"Merging clusters {clusters_from} into '{cluster_to}'...")

        # Count cells before merging
        counts_from = {}
        for cluster in clusters_from:
            counts_from[cluster] = (adata.obs[new_resolution] == cluster).sum()

        count_to = (adata.obs[new_resolution] == cluster_to).sum()
        
        print(f"Cells in '{cluster_to}' before merging: {count_to}")
        for cluster, count in counts_from.items():
            print(f"  - Cluster '{cluster}': {count} cells")

        # Update the labels in the new resolution column
        replacement_dict = {}
        for cluster in clusters_from:
            replacement_dict[cluster] = cluster_to

        adata.obs[new_resolution] = adata.obs[new_resolution].replace(replacement_dict)

        # Count cells after merging
        new_count_to = (adata.obs[new_resolution] == cluster_to).sum()
        print(f"Cells in '{cluster_to}' after merging: {new_count_to}")

    return adata

def print_cell_counts(adata, resolution):
    """
    Print the number of cells in each cluster for a given resolution.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    resolution (str): The resolution column in `adata.obs` to analyze.

    Returns:
    None
    """

    print(f"Cell counts for resolution '{resolution}':\n")
    cell_counts = adata.obs[resolution].value_counts()

    for cluster, count in cell_counts.items():
        print(f"Cluster '{cluster}': {count} cells")

def save_adata(adata, file_path):
    """
    Save the AnnData object to an .h5ad file with gzip compression.

    Parameters:
    adata (AnnData): The AnnData object to save.
    file_path (str): The file path to save the .h5ad file.

    Returns:
    None
    """
    t0 = time.time()  # Start timing
    print(f"Saving AnnData object to '{file_path}'...")
    
    # Save the file
    adata.write_h5ad(file_path, compression='gzip')
    
    print(f"File saved in '{file_path}'. Time elapsed: {time.time() - t0:.1f} seconds.")


# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")

    merge_groups = [
    (['MeV.4.4', 'MeV.1.4.0'], 'MeV.4.4'),
    (['MeV.1.4.12', 'MeV.NA'], 'MeV.NA')  
]

    adata = merge_clusters(adata, 'leiden_fusion_old2', merge_groups, 'leiden_fusion')

    print_cell_counts(adata, 'leiden_fusion')

    output_file = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    save_adata(adata, output_file)