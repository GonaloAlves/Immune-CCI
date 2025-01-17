# Import necessary packages
import os
import scanpy as sc
import time

# Step 1: Load the dataset
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

# Step 2: Perform DGE analysis
def dge_data(adata, groupby, uns_key):
    """
    Perform DGE analysis and store the ranked gene groups in the `uns` attribute.

    Parameters:
    adata (AnnData): The AnnData object.
    groupby (str): The column in `adata.obs` to group cells by for DGE analysis.
    uns_key (str): The key name to store ranked gene data in `adata.uns`.

    Returns:
    AnnData: The updated AnnData object with ranked gene data stored in `uns`.
    """
    print(f"Performing DGE analysis using groupby='{groupby}'...")
    sc.tl.rank_genes_groups(adata, groupby, method='wilcoxon', use_raw=False, pts=True)
    
    # Save the results in `uns` under the specified key
    adata.uns[uns_key] = adata.uns['rank_genes_groups']

    if 'rank_genes_groups_leiden_mako' in adata.uns:
        print(f"Ranked genes stored in `uns` with key '{uns_key}'.")
    
    return adata

# Step 3: Save the AnnData object
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
    file_path = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    adata = load_data(file_path)

    # Perform DGE analysis
    adata = dge_data(adata, 'leiden_mako', 'rank_genes_groups_leiden_mako')

    # Save the updated AnnData object
    output_file = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    save_adata(adata, output_file)
