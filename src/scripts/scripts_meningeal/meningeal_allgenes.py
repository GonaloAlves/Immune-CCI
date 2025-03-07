# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Load the dataset
def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    return sc.read_h5ad(file_path)




def export_top_genes_to_txt(top_genes_cluster, threshold, output_dir="excels/meningeal/updates"):
    """
    Export the top genes from each cluster to a single text file.

    Parameters:
    top_genes_cluster (dict): Dictionary containing top genes for each cluster.
    threshold (float): The threshold value used for filtering, included in the filename.
    output_dir (str): Directory to save the text file.

    Returns:
    None
    """
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Define the output file path
    output_file = os.path.join(output_dir, f"top_genes_clusters_{threshold}.txt")

    print(f"\nExporting top genes to text file: {output_file}")

    with open(output_file, 'w') as f:
        for cluster, df in top_genes_cluster.items():
            f.write(f"Cluster: {cluster}\n")  # Write cluster name
            f.write(f"Top Genes:\n")
            for gene in df.index.tolist():
                f.write(f"  - {gene}\n")  # Write each gene name
            f.write("\n" + "-"*40 + "\n")  # Separator between clusters

    print(f"Text file saved: {output_file}")