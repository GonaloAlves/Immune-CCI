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
    print("dge_immune.py")
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

    if 'rank_genes_groups_leiden_fusion_old1' in adata.uns:
        print(f"Ranked genes stored in `uns` with key '{uns_key}'.")
    
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


def dendogram_sc(adata):
    """
    Compute and save the dendrogram for a specific group and store it in a specified key.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    dendrogram_key (str): The key to save the dendrogram in `adata.uns`.

    Returns:
    None
    """

    # Compute the dendrogram
    print(f"Computing dendrogram for leiden_fusion...")
    sc.tl.dendrogram(
        adata,
        groupby='leiden_fusion_old1',
        use_rep='X_pca',
        cor_method='spearman',
        linkage_method='ward',
        use_raw=False
    )

def remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    mask_NA = adata.obs['leiden_fusion_old1'] != 'Imm.NA' #creates mask for remove NA cells
    #print(mask_NA)    
    adata2 = adata[mask_NA] #apply mask

    # print(adata2.obs['leiden_fusion'].cat.categories.to_list())
    # adata2.obs['leiden_fusion'] = adata2.obs['leiden_fusion'].cat.remove_unused_categories()
    # print(adata2.obs['leiden_fusion'].cat.categories.to_list())
    # print(adata.obs['leiden_fusion'])

    return adata2

# Main execution block
if __name__ == "__main__":
    # Load data
    file_path = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad"
    adata = load_data(file_path)

    # Perform DGE analysis
    adata = dge_data(adata, 'leiden_fusion_old1', 'rank_genes_groups_leiden_fusion_old1')

    #filtered_adata = remove_NA_cat(adata)

    # Preform Dendrogram
    dendogram_sc(adata)
        
    # adata = drop_mako(adata)
    # print(adata)
    print(adata.obs['leiden_fusion_old1'].cat.categories.to_list())

    # Save the updated AnnData object
    output_file = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad"
    save_adata(adata, output_file)
