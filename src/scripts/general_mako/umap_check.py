
import os
import shutil
import scanpy as sc
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

def load_data(file_path):
    """
    Load the data from a .h5ad file.
    
    Parameters:
    file_path (str): Path to the .h5ad file.

    Returns:
    AnnData: The loaded AnnData object.
    """
    print(f"Loading h5ad file: {file_path}...")

    return sc.read_h5ad(file_path)

def load_obsm_from_excel(adata, excel_path):
    """
    Load obsm matrices from an Excel file and add them to an AnnData object.

    Parameters:
    adata (AnnData): The AnnData object to update.
    excel_path (str): Path to the Excel file containing obsm data.
    """
    print(f"Loading obsm data from: {excel_path}...")
    obsm_data = pd.read_excel(excel_path, sheet_name=None, index_col=0)
    
    for key, df in obsm_data.items():
        df = df.reindex(adata.obs.index)
        adata.obsm[key] = df.values
    
    adata.obsm["X_umap"] = adata.obsm["X_umap_leiden_fusion"]
    
    print("obsm data successfully loaded.")

def load_obsp_from_excel(adata, excel_path):
    """
    Load obsp matrices from an Excel file and add them to an AnnData object.

    Parameters:
    adata (AnnData): The AnnData object to update.
    excel_path (str): Path to the Excel file containing obsp data.
    """
    print(f"Loading obsp data from: {excel_path}...")
    obsp_data = pd.read_excel(excel_path, sheet_name=None, index_col=0)
    
    for key, df in obsp_data.items():
        df = df.reindex(index=adata.obs.index, columns=adata.obs.index)
        adata.obsp[key] = csr_matrix(df.values)
    
    print("obsp data successfully loaded.")


def set_umap_key(adata):
    """
    Set the correct UMAP key and update the neighbors for Leiden fusion clustering.
    
    Parameters:
    - adata (AnnData): Annotated data matrix.
    
    Returns:
    - None
    """
    print("Updating UMAP and neighbors for Leiden Fusion...")
    
    # Set the correct UMAP embedding
    if "X_umap_leiden_fusion" in adata.obsm:
        adata.obsm["X_umap"] = adata.obsm["X_umap_leiden_fusion"]
        print("✔ UMAP key updated to 'X_umap_leiden_fusion'")
    else:
        print("⚠ WARNING: 'X_umap_leiden_fusion' not found in adata.obsm")

    # Set the correct PCA embedding
    if "X_pca_leiden_fusion" in adata.obsm:
        adata.obsm["X_pca"] = adata.obsm["X_pca_leiden_fusion"]
        print("✔ PCA key updated to 'X_pca_leiden_fusion'")
    else:
        print("⚠ WARNING: 'X_pca_leiden_fusion' not found in adata.obsm")
    
    # Update the neighbors key if needed
    if "neighbors_15_connectivities_leiden_fusion" in adata.obsp:
        adata.obsp["neighbors_15_connectivities"] = adata.obsp["neighbors_15_connectivities_leiden_fusion"]
        adata.obsp["neighbors_15_distance"] = adata.obsp["neighbors_15_distances_leiden_fusion"]
        print("✔ Neighbors updated for Leiden Fusion")
    else:
        print("⚠ WARNING: Leiden fusion neighbors not found in adata.obsp")

def compare_umap(adata_immune, adata_meningeal ,output_dir='/home/makowlg/Documents/Immune-CCI/src/reso/resos'):
    """
    Compare and plot UMAP embeddings from two datasets.
    
    Parameters:
    adata_immune (AnnData): AnnData object for the Immune dataset.
    adata_meningeal (AnnData): AnnData object for the Meningeal Vascular dataset.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    sc.pl.umap(adata_immune, color='leiden_fusion', ax=axes[0], title='Immune Dataset')
    sc.pl.umap(adata_meningeal, color='leiden_fusion', ax=axes[1], title='Meningeal Vascular Dataset')
    
    output_path = os.path.join(output_dir, 'umap_comparison.png')
    plt.savefig(output_path, bbox_inches='tight')
    print(f"UMAP comparison saved to {output_path}")

def umap_reso_cluster(adata, resolution_name, keyname, d, output_dir="reso/resos"):
    """
    Plot UMAP for a specific resolution with cluster numbers displayed at the centroid of each cluster.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    resolution_name (str): The resolution name to be used for coloring the UMAP plot (e.g., 'leiden_fusion').

    Returns:
    None
    """

    # Plot UMAP for the specified resolution
    ax = sc.pl.umap(adata, color=resolution_name ,title=f"UMAP - {resolution_name}", return_fig=True)
    ax_ondata = sc.pl.umap(adata, color=resolution_name, title=f"UMAP - {resolution_name}",legend_loc = 'on data',legend_fontweight = 'bold', legend_fontsize = 12, legend_fontoutline = 1,return_fig=True)

    # Save the UMAP plot as an image (optional)
    output_path= os.path.join(output_dir, f"umap_{d}_{keyname}_{resolution_name}_n.png") 
    output_dir_leg = os.path.join(output_dir, f"umap_{d}_{keyname}_{resolution_name}_l.png")
    ax.figure.savefig(output_path, bbox_inches="tight")
    ax_ondata.figure.savefig(output_dir_leg, bbox_inches="tight")
    # print(f"UMAP plot saved as {output_path}")
    plt.close()  # Close the plot to avoid overlap

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
    

if __name__ == "__main__":
    datasets = ["Immune", "Meningeal_Vascular"]
    adatas = {}

    for d in datasets:
        file_path = f"/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_{d}_raw_norm_ranked_copy_copy.h5ad"
        excel_path_obsm = f"/home/makowlg/Documents/Immune-CCI/src/reso/umap_cords/{d}_obsm_data.xlsx"
        excel_path_obsp = f"/home/makowlg/Documents/Immune-CCI/src/reso/umap_cords/{d}_obsp_data.xlsx"
        
        adata = load_data(file_path)
        umap_reso_cluster(adata,  "leiden_fusion", "old", d)

        print(adata)
        load_obsm_from_excel(adata, excel_path_obsm)
        load_obsp_from_excel(adata, excel_path_obsp)
        set_umap_key(adata)
        adatas[d] = adata
        print(adata)
    
        # umap_reso_cluster(adata, "leiden_fusion", "new", d)

        # Save the updated AnnData object
        output_file = f"/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_{d}_raw_norm_ranked_copy_copy.h5ad"
        save_adata(adata, output_file)



