# Import packages
import os
import scanpy as sc
import time


# Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    print("Loading h5ad file...")
    return sc.read_h5ad(file_path)


def rename_clusters(adata, rename_pairs, resolution="leiden_fusion"):
    """
    Rename multiple clusters in the provided resolution column.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    rename_pairs (list of tuples): A list of (current_cluster_name, new_cluster_name) pairs.
    resolution (str): The resolution column in `adata.obs` to modify (default: "leiden_fusion").

    Returns:
    AnnData: The modified AnnData object with renamed clusters.
    """

    # Rename clusters
    for current_name, new_name in rename_pairs:
        if current_name not in adata.obs[resolution].unique():
            print(f"Warning: Cluster '{current_name}' not found in resolution '{resolution}'. Skipping...")
            continue

        print(f"Renaming cluster '{current_name}' to '{new_name}' in resolution '{resolution}'...")
        adata.obs[resolution] = adata.obs[resolution].replace({current_name: new_name})

        # Confirm the rename
        print(f"Cluster '{current_name}' has been renamed to '{new_name}'.")
        print(adata.obs[resolution].value_counts())

    return adata


# Main execution block
if __name__ == "__main__":
    # Load data
    file_path = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    adata = load_data(file_path)

    # List of clusters to rename
    old_clusters = [
        ('Imm.0.8.6', 'Imm.Interferon'),
        ('Imm.0.8.0', 'Imm.M0_like1'),
        ('Imm.0.8.1', 'Imm.M0_like2'),
        ('Imm.0.8.2', 'Imm.M0_like3'),
        ('Imm.1.2.15', 'Imm.Proliferation'),
        ('Imm.1.2.4', 'Imm.PVM.0'),
        ('Imm.1.2.12', 'Imm.DAM.0'),
        ('Imm.1.2.5', 'Imm.DAM.1')
    ]
    
    rename_pairs = [
        ('MeV.1.4.1', 'MeV.Vascular.0'),
        ('MeV.4.21', 'MeV.Vascular.1'),
        ('MeV.1.4.5', 'MeV.Vascular.2'),
        ('MeV.4.31', 'MeV.SMC.0'),
        ('MeV.4.1', 'MeV.Pericytes.0'),
        ('MeV.3.30', 'MeV.Proliferative_Fibr.0')
    ]

    # Rename the clusters
    adata = rename_clusters(adata, rename_pairs)

    print(adata.obs['leiden_fusion'])

    # Save the modified AnnData object
    output_path = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    print(f"Saving modified AnnData to '{output_path}'...")
    adata.write_h5ad(output_path, compression="gzip")
    print("Save complete.")
