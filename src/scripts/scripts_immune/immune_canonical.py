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
    
    mask_NA = adata.obs['leiden_mako'] != 'Imm.NA' #creates mask for remove NA cells 
    adata2 = adata[mask_NA] #apply mask
    return adata2

# Load canonical genes
def load_canonical_from_dir(directory):
    """
    Load gene lists from all .txt files in a directory.

    Parameters:
    directory (str): Path to the directory containing .txt files.

    Returns:
    dict: Dictionary of gene lists with filenames (without extensions) as keys.
    """
    print(f"Loading canonical gene lists from: {directory}")
    gene_files = [f for f in os.listdir(directory) if f.endswith('.txt')]

    gene_dict = {}
    for gene_file in gene_files:
        file_path = os.path.join(directory, gene_file)
        gene_name = os.path.splitext(gene_file)[0]  # Use the filename without the extension
        with open(file_path, 'r') as f:
            gene_dict[gene_name] = f.read().splitlines()

    print(f"Loaded gene lists: {list(gene_dict.keys())}")
    return gene_dict

# Step 7: Visualize the DotPlots of the DGE's
def create_dotplot(adata, top_genes_names, output_dir="canonical_immune"):
    """
    Create and save a dotplot of the top genes per cluster, ordered alphabetically by gene group.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    top_genes_names (dict): Dictionary of top gene names per cluster.
    output_dir (str): Directory to save the dotplot.

    Returns:
    None
    """
    # Create the directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Sort the gene groups (keys) alphabetically
    top_genes_names = {key: top_genes_names[key] for key in sorted(top_genes_names.keys())}

    # Generate the dotplot
    print(f"\nGenerating a dotplot")
    dotplot = sc.pl.dotplot(
        adata,
        var_names=top_genes_names,
        groupby='leiden_mako',
        cmap='Greys',
        vmin=0,
        vmax=1,
        colorbar_title='Mean expression',
        use_raw=False,
        standard_scale='var',
        #dendrogram = False,
        dendrogram='dendrogram_leiden_mako',
        return_fig=True
    )

    # Save the plot
    output_path = os.path.join(output_dir, "dotplot_immune_canonical_ordered_dani_dendro.png")
    dotplot.savefig(output_path, bbox_inches="tight")
    plt.close()  # Close the current figure to avoid overlap


def dendogram_sc(adata):
    
    """
    
    """
    # Compute the dendrogram
    print(f"Computing dendrogram for leiden_mako...")
    sc.tl.dendrogram(
        adata,
        groupby='leiden_mako',
        use_rep= 'X_pca',
        cor_method= 'spearman',
        linkage_method='ward',
        use_raw=False
    )



# Main execution block
if __name__ == "__main__":
    # Load data
    adata = load_data("/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad")

    filtered_adata = remove_NA_cat(adata)

    # Load canonical gene lists from a directory
    canonical_genes_dir = "/home/makowlg/Documents/Immune-CCI/src/canonical_txt/Immune"
    genes = load_canonical_from_dir(canonical_genes_dir)

    # Create dendogram ot the top genes
    dendogram_sc(filtered_adata)

    # Create dotplot of the top genes
    create_dotplot(filtered_adata, genes)

    