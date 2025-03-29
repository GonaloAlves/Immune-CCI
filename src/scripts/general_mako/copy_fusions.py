import os
import scanpy as sc
import pandas as pd
import time

def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    print("copy_fusions.py")
    print("Loading h5ad file")
    return sc.read_h5ad(file_path)


def copy_leidens(adata, old_reso1, new_reso1, fusion):
    """
    Copy and organize leiden resolutions in the AnnData object.

    Parameters:
    adata (AnnData): The AnnData object containing the resolutions.
    old_reso1 (str): The current name of the first resolution to duplicate.
    new_reso1 (str): The new name for the first resolution copy.
    old_reso2 (str): The current name of the second resolution to duplicate.
    new_reso2 (str): The new name for the second resolution copy.
    fusion (str): The target resolution name for the second resolution.

    Returns:
    AnnData: The modified AnnData object with the new resolution columns.
    DataFrame: A DataFrame containing the unique clusters for comparison.
    """
    print("Duplicating resolutions into the correct forms...")
    print(f"{old_reso1} ---> {new_reso1}")
    #print(f"{old_reso2} ---> {new_reso2}")
    print(f"{old_reso1} ---> {fusion}")
    
    # Create the new resolution columns by copying existing ones
    adata.obs[old_reso1] = adata.obs[new_reso1]  # Copy old_reso1 to new_reso1
    #adata.obs[new_reso2] = adaold_reso1ta.obs[old_reso2]  # Copy old_reso2 to new_reso2
    adata.obs[fusion] = adata.obs[new_reso1]     # Copy fusion to old_reso1
    
    # # Remove the 'leiden_mako' column (or `old_reso2`) from adata.obs
    # adata.obs.drop(columns=[old_reso2], inplace=True)
    # print(f"Removed '{old_reso2}' from adata.obs.")
    
    # Extract unique clusters for each remaining resolution
    unique_clusters = {
        new_reso1: adata.obs[new_reso1].unique().tolist(),
        #new_reso2: adata.obs[new_reso2].unique().tolist(),
        fusion: adata.obs[fusion].unique().tolist()
    }

    print("Unique clusters in each resolution:")
    for reso, clusters in unique_clusters.items():
        print(f"{reso}: {clusters}")

    # Convert unique clusters to a DataFrame for saving
    txt = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in unique_clusters.items()]))

    return adata, txt


def save_txt(txt, file_path="unique_clusters_comparison.xlsx"):
    """
    Save the unique clusters to an Excel file for comparison.

    Parameters:
    txt (DataFrame): DataFrame containing the unique clusters to save.
    file_path (str): The file path to save the Excel file.

    Returns:
    None
    """
    print(f"Saving unique clusters to '{file_path}'...")
    txt.to_excel(file_path, index=False)  # Save without index
    print(f"File saved as '{file_path}'.")


# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad")
    print("leiden_fusion")
    print(adata.obs['leiden_fusion'].cat.categories.to_list())

    # Copy and organize leiden resolutions
    adata, txt = copy_leidens(
        adata, 
        old_reso1='leiden_fusion_old1', 
        new_reso1='leiden_fusion', 
        # old_reso2='leiden_mako', 
        # new_reso2='leiden_fusion_old2', 
        fusion='leiden_fusion'
    )
    print(adata)
    print(adata.obs['leiden_fusion'].cat.categories.to_list())
    #print(adata.obs['annotation'].cat.categories.to_list())
    print(adata.obs['leiden_fusion_old1'].cat.categories.to_list())
    # Save the unique clusters to an Excel file
    save_txt(txt, file_path="unique_clusters_comparison.xlsx")

    # Save the updated AnnData object (optional)
    output_path = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Meningeal_Vascular_raw_norm_ranked_copy_copy.h5ad"
    print(f"Saving updated AnnData to {output_path}...")
    adata.write_h5ad(output_path, compression='gzip')
    print("Save complete.")
