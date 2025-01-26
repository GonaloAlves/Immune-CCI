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

# Step 2: Backup DGE data in `uns`
def backup_dge_data(adata, old_key, backup_key):
    """
    Backup the existing DGE data from `old_key` to `backup_key` in `uns`.

    Parameters:
    adata (AnnData): The AnnData object.
    old_key (str): The key in `adata.uns` to back up.
    backup_key (str): The key to store the backup in `adata.uns`.

    Returns:
    AnnData: The updated AnnData object with the backed-up data.
    """
    if old_key in adata.uns:
        print(f"Backing up '{old_key}' to '{backup_key}' in `uns`...")
        adata.uns[backup_key] = adata.uns[old_key]
        print(f"Backup complete: '{backup_key}' now contains the previous DGE data.")
    else:
        print(f"Warning: Key '{old_key}' not found in `uns`. No backup created.")
    return adata


def drop_mako(adata, dge_key='rank_genes_groups_leiden_mako'):
    """
    Remove the 'leiden_mako' column and its corresponding DGE key from AnnData.

    Parameters:
    adata (AnnData): The AnnData object.
    old_reso2 (str): The column name in `adata.obs` to remove (e.g., 'leiden_mako').
    dge_key (str): The key in `adata.uns` to remove (default: 'rank_genes_groups_leiden_mako').

    Returns:
    AnnData: The updated AnnData object with the specified column and key removed.
    """

    # Remove the 'rank_genes_groups_leiden_mako' key from adata.uns
    if dge_key in adata.uns:
        del adata.uns[dge_key]
        print(f"Removed '{dge_key}' from adata.uns.")
    else:
        print(f"'{dge_key}' not found in adata.uns. No action taken.")

    return adata


# Step 3: Perform DGE analysis
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

    print(f"Ranked genes stored in `uns` with key '{uns_key}'.")
    return adata

# Step 4: Save the AnnData object
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

    # Backup the existing DGE data
    adata = backup_dge_data(adata, 
                             old_key='rank_genes_groups_leiden_fusion', 
                             backup_key='rank_genes_groups_leiden_fusion_old1')

    # Perform new DGE analysis
    adata = dge_data(adata, 'leiden_fusion', 'rank_genes_groups_leiden_fusion')

    #data = drop_mako(adata)


    print(adata)
    # Save the updated AnnData object
    output_file = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    save_adata(adata, output_file)
