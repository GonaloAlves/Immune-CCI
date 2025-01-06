# Import necessary packages
import os
import shutil
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# This script will display the gene expression values of the canonical genes from Immune dataset


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


def remove_NA_cat(adata: sc.AnnData):
    
    print("Removing NA cells category")
    
    mask_NA = adata.obs['leiden_fusion'] != 'Imm.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2

# Load canonical genes
def load_canonical(genes_txt):

    genes = open(genes_txt).read().split()
    print(genes)

    top_genes_names = {'Immune': genes }
    print (top_genes_names)

    return top_genes_names

# Step 7: Visualize the DotPlots of the DGE's
def create_dotplot(adata, top_genes_names, output_dir="canonical_immune"):
    """
    Create and save a dotplot of the top genes per cluster.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    top_genes_names (dict): Dictionary of top gene names per cluster.
    output_dir (str): Directory to save the dotplot.

    Returns:
    None
    """
    # Create the directory if it doesn't exist
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Generate the dotplot
    print(f"\nGenerating a dotplot")

    dotplot = sc.pl.dotplot(
        adata,
        var_names=top_genes_names,
        groupby='leiden_fusion',
        cmap='bwr',
        vmin=-4,
        vmax=4,
        colorbar_title='log fold change',
        use_raw=False,
        dendrogram=False,
        #dendrogram='dendrogram_leiden_fusion',
        return_fig=True
    )

    output_path = os.path.join(output_dir, "dotplot_immune_canonical.png")
    dotplot.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap



def dendogram_sc(adata):
    
    """
    
    """
    # Compute the dendrogram
    print(f"Computing dendrogram for leiden_fusion...")
    sc.tl.dendrogram(
        adata,
        groupby='leiden_fusion',
        use_rep= 'X_pca',
        cor_method= 'spearman',
        linkage_method='ward',
        use_raw=False
    )



# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad")

    filtered_adata = remove_NA_cat(adata)

    genes = load_canonical("/home/makowlg/Documents/Immune-CCI/src/Immune_canonical_genes.txt")

    # Create dendogram ot the top genes
    dendogram_sc(filtered_adata)

    # Create dotplot of the top genes
    create_dotplot(filtered_adata, genes)

    