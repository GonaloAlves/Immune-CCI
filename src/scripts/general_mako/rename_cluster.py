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
    print("rename_cluster.py")
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
        ('Imm.0.8.6', 'Imm.Interferon.0'),
        ('Imm.0.8.0', 'Imm.M0Like1'),
        ('Imm.0.8.1', 'Imm.M0Like2'),
        ('Imm.0.8.2', 'Imm.M0Like3'),
        ('Imm.1.2.15', 'Imm.Proliferative.0'),
        ('Imm.1.2.4', 'Imm.PVM.0'),
        ('Imm.1.2.12', 'Imm.DAM.0'),
        ('Imm.1.2.5', 'Imm.DAM.1'),
        ('Imm.1.2.13', "Imm.MHCII.0")
    ]
    
    old_meningeal = [
        ('MeV.1.4.1', 'MeV.Endothelial.0'),
        ('MeV.4.21', 'MeV.Endothelial.3'),
        ('MeV.1.4.5', 'MeV.Endothelial.2'),
        ('MeV.1.4.15', 'MeV.Endothelial.1'),
        ('MeV.1.4.20', 'MeV.EndoUnknow.4'),
        ('MeV.4.31', 'MeV.SMC.0'),
        ('MeV.4.1', 'MeV.Pericytes.0'),
        ('MeV.3.30', 'MeV.FibProlif.0'),
        ('MeV.4.30', 'MeV.FibCollagen.2'),
        ('MeV.1.4.8', 'MeV.LowQuality.0'),
        ('MeV.1.4.21', 'MeV.ImmuneDoublets.0'),
        ('MeV.4.26', 'MeV.Epithelial.0'),
        ('MeV.4.4', 'MeV.VLMC.0'),
        ('MeV.4.12', 'MeV.VLMC.1'),
        ('MeV.3.17', 'MeV.FibLaminin.0'),
        ('MeV.2.1', 'MeV.Fib.0'),
        ('MeV.2.8', 'MeV.FibCollagen.0'),
        ('MeV.1.4.2', 'MeV.FibUnknown.6'),
        ('MeV.1.4.11', 'MeV.Fib.5'),
        ('MeV.1.4.4', 'MeV.FibCollagen.3'),
        ('MeV.1.4.7', 'MeV.Fib.1'),
        ('MeV.1.4.6', 'MeV.Fib.2'),
        ('MeV.1.4.12', 'MeV.FibCollagen.1'),
        ('MeV.4.34', 'MeV.Fib.3'),
        ('MeV.1.4.13', 'MeV.Fib.4')
    ]

    rename_pairs= [
        ('MeV.1.4.1', 'MeV.Endothelial.0'),
        ('MeV.4.21', 'MeV.Endothelial.3'),
        ('MeV.1.4.5', 'MeV.Endothelial.2'),
        ('MeV.1.4.15', 'MeV.Endothelial.1'),
        ('MeV.1.4.20', 'MeV.EndoUnknow.4'),
        ('MeV.4.31', 'MeV.SMC.0'),
        ('MeV.4.1', 'MeV.Pericytes.0'),
        ('MeV.3.30', 'MeV.FibProlif.0'),
        ('MeV.4.30', 'MeV.FibCollagen.2'),
        ('MeV.1.4.8', 'MeV.LowQuality.0'),
        ('MeV.1.4.21', 'MeV.ImmuneDoublets.0'),
        ('MeV.4.26', 'MeV.Epithelial.0'),
        ('MeV.4.4', 'MeV.VLMC.0'),
        ('MeV.4.12', 'MeV.VLMC.1'),
        ('MeV.3.17', 'MeV.FibLaminin.0'),
        ('MeV.2.1', 'MeV.Fib.0'),
        ('MeV.2.8', 'MeV.FibCollagen.0'),
        ('MeV.1.4.2', 'MeV.FibUnknown.6'),
        ('MeV.1.4.11', 'MeV.Fib.5'),
        ('MeV.1.4.4', 'MeV.FibCollagen.3'),
        ('MeV.1.4.7', 'MeV.Fib.1'),
        ('MeV.1.4.6', 'MeV.Fib.2'),
        ('MeV.1.4.12', 'MeV.FibCollagen.1'),
        ('MeV.4.34', 'MeV.Fib.3'),
        ('MeV.1.4.13', 'MeV.Fib.4')
    ]

    # Rename the clusters
    adata = rename_clusters(adata, rename_pairs)

    print("----")
    print(adata.obs['leiden_fusion'].cat.categories.to_list())
    print("----")

    # Save the modified AnnData object
    output_path = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    print(f"Saving modified AnnData to '{output_path}'...")
    adata.write_h5ad(output_path, compression="gzip")
    print("Save complete.")
