# Perform CellPhoneDB analysis - LeonorSaude10x
#
# Follwing these tutorials
#
# Daniel Ribeiro, Gonçalo Alves 2025
import gc
from typing import Optional  # Import Optional
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('AGG')   # Use non-interactive backend for plotting
import matplotlib.pyplot as plt
import ktplotspy as kpy
import seaborn as sns



statistical_analysis = True
deg_analysis = True

checkpoint_dir = "/home/makowlg/Documents/Immune-CCI/h5ad_files"
cellphonedb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb"

def get_cell_types(cci: str,
                   major: bool      # Only major cell types, no clusters
                   
                   ) -> tuple[str, ...]:
    cell_types = cci.split('|')
    if major:
        cell_types[0] = cell_types[0].split('.')[0]
        cell_types[1] = cell_types[1].split('.')[0]
    
    return tuple(cell_types)

def plot_heatmaps(adata: sc.AnnData, log1p: bool = False) -> None:
    """
    Plots a heatmap using ktplotspy for a given AnnData object.

    Parameters
    ----------
    `adata` : sc.AnnData
        AnnData object used for plotting.
        
    `log1p` : bool, optional
        Whether to log-transform the number of significant interactions for better visualization.
    """

    # Load p-values file (assuming a fixed filename)
    dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged.txt"
    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)

    # Generate heatmap
    clusterg = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        log1p_transform=log1p,
        figsize=(60, 60),
        linewidths=0.2
    )

    # Title customization
    suptitle = "Number of significant interactions (log1p transformed)" if log1p \
        else "Number of significant interactions"
    clusterg.figure.suptitle(suptitle, y=0.85, size=50)
    clusterg.ax_cbar.set_position((1.02, 0.0, .02, .3))
    
    # Adjust heatmap
    ax = clusterg.ax_heatmap
    ax.grid(False)
    
    # Save plot
    output_path = f"{cellphonedb_dir}/significant_interactions_final_merged.pdf"
    print(f"Saving heatmap to: {output_path}")
    clusterg.savefig(output_path, bbox_inches="tight")
    plt.close()

def plot_heatmaps_major_cells(obs_key: Optional[str] = None,
                              category: Optional[str] = None,
                              log1p: bool = False) -> None:
    """
    Parameters
    ----------
    `obs_key` : str
        Key to retreive available categories. Not fully implemented, except to evaluate whether \
            somenone asks for a category without providing a obs_key.
    
    `category` : str
        Category name to be found in .txt file name.
    
    `log1p` : bool
        Wether to log the number of significant categories. May improve heatmap vizualization.
    """
    
    if obs_key is not None and category is None:
        raise ValueError("If obs_key is not None then category cannot be None!")
    elif obs_key is not None:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{category.replace('.', '_')}.txt"
    else:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged.txt"
        
    pvalues = pd.read_csv(dest, sep='\t', index_col=0, dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
    pvalues.drop(columns=["interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted",
                          "receptor_a", "receptor_b", "annotation_strategy", "is_integrin"], inplace=True)
    
    # Build a source/target matrix - row = source, col = target
    # Get all possible sources and targets
    sources = []
    targets = []
    for col in pvalues.columns:
        cci = get_cell_types(col, major=True)
        sources.append(cci[0])
        targets.append(cci[1])
    sources = list(set(sources))
    sources.sort()
    targets = list(set(targets))
    targets.sort()
    mtx_source_target = pd.DataFrame(np.zeros((len(sources), len(targets))), index=sources, columns=targets, dtype='int64')
    mtx_source_target.index.name = "source"
    mtx_source_target.columns.name = "target"
    
    # For every source/target significant interaction add a +1 weight
    # Get all ccis for each src/tgt combination
    for col in pvalues.columns:
        mask = pvalues.loc[:, col] < 0.05
        n = mask.sum()
        src, tgt = get_cell_types(col, major=True)
        mtx_source_target.loc[src, tgt] += n
    
    if log1p:
        mtx_source_target = np.log1p(mtx_source_target)
    
    # Plot heatmap
    cm = sns.color_palette("rocket_r", as_cmap=True)
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    ax = sns.heatmap(mtx_source_target,
                     cmap=cm,
                     linewidths=0.5,
                     square=True,
                     cbar=True,
                     cbar_kws={'pad': 0.1, 'shrink': 0.7, 'label': "log1p significant interactions"})
    ax.grid(False)
    suptitle = "Number of significant interactions for major cell types (log1p transformed)" if log1p \
        else "Number of significant interactions for major cell types"
    plt.title(suptitle)
    if obs_key is not None:
        dest = cellphonedb_dir + f"/significant_interactions_final_merged_major_{category.replace('.', '_')}.pdf"
    else:
        dest = cellphonedb_dir + "/significant_interactions_final_merged_major.pdf"
    print(dest)
    plt.savefig(dest, bbox_inches='tight')
    plt.close(fig)
    
    # Save mtx_source_target matrix
    dest = f"mtx_source_target_{category.replace('.', '_')}" if obs_key is not None else \
        "mtx_source_target"
    dest = f"{cellphonedb_dir}/{dest}_log1p.txt" if log1p else \
        f"{cellphonedb_dir}/{dest}.txt"
    mtx_source_target.to_csv(dest, sep='\t')




def plot_lineage_vs_other_interactions(adata: sc.AnnData, lineage_prefix: str) -> None:

    # Check if lineage exists
    cell_types = adata.obs["leiden_merge"].cat.categories.str.startswith(lineage_prefix)
    if not cell_types.any():
        raise ValueError(f"Lineage prefix '{lineage_prefix}' not found!")

    # Load interaction data (ignoring injury conditions)
    dest_pvalues = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged.txt"
    dest_means = f"{cellphonedb_dir}/statistical_analysis_means_final_merged.txt"
    
    pvalues = pd.read_csv(dest_pvalues, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    means = pd.read_csv(dest_means, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})

    # Define cell groups
    all_cell_types = adata.obs["leiden_merge"].cat.categories.to_list()
    cell_types1 = [cell for cell in all_cell_types if cell.startswith(lineage_prefix)]
    cell_types2 = [cell for cell in all_cell_types if not cell.startswith(lineage_prefix)]
    
    # Ensure lists are not empty
    if not cell_types1:
        raise ValueError(f"No cell types found with prefix '{lineage_prefix}'!")
    if not cell_types2:
        raise ValueError("No other cell types found!")

    

    # Join cell type names into strings
    cell_types1 = "|".join(cell_types1)
    cell_types2 = "|".join(cell_types2)

    # Create dot plot
    dotplot = kpy.plot_cpdb(
        adata=adata,
        cell_type1=cell_types1,
        cell_type2=cell_types2,
        means=means,
        pvals=pvalues,
        celltype_key="leiden_merge",
        figsize=(40, 30),
        title=f"Interactions between '{lineage_prefix}' and other cells",
        max_size=3,
        highlight_size=1,
        standard_scale=True
    )

    # Save plot
    dest_plot = f"{cellphonedb_dir}/{lineage_prefix}_interactions_final_merged.pdf"
    dotplot.save(dest_plot, dpi=300, limitsize=False, bbox_inches="tight")
    plt.close()

    print(f"Saved plot: {dest_plot}")

def chord_diagram(adata: sc.AnnData) -> None:

    # Load interaction data (ignoring injury conditions)
    dest_pvalues = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged.txt"
    dest_means = f"{cellphonedb_dir}/statistical_analysis_means_final_merged.txt"

    pvalues = pd.read_csv(dest_pvalues, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    means = pd.read_csv(dest_means, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})

    all_cell_types = adata.obs["leiden_merge"].cat.categories.to_list()

    # Use all cell types for both groups to compare everything against everything
    cell_types1 = "|".join(all_cell_types)
    cell_types2 = "|".join(all_cell_types)

    deconvoluted_file = f"{cellphonedb_dir}/statistical_analysis_deconvoluted_final_merged.txt"
    deconvoluted_data = pd.read_csv(deconvoluted_file, sep="\t")

    #Debugging

    # Check cell types in adata
    adata_cell_types = set(adata.obs["leiden_merge"].unique())

    # Check cell types in the deconvoluted data
    deconvoluted_cell_types = set(deconvoluted_data["interacting_pair"].unique())

    # Find mismatches
    only_in_adata = adata_cell_types - deconvoluted_cell_types
    only_in_deconvoluted = deconvoluted_cell_types - adata_cell_types

    print(f"Cell types only in adata: {only_in_adata}")
    print(f"Cell types only in deconvoluted data: {only_in_deconvoluted}")



    chord = kpy.plot_cpdb_chord(
        adata=adata,
        cell_type1=cell_types1,
        cell_type2=cell_types2,
        means=means,
        pvals=pvalues,
        deconvoluted= deconvoluted_data,
        celltype_key="leiden_merge",
        figsize=(40, 30),
        title="Interactions between all cells",
        max_size=3,
        highlight_size=1,
        standard_scale=True
    )

    
    # Save plot
    dest_plot = f"{cellphonedb_dir}/chord_all_interactions_final_merged.pdf"

    # If `chord` is a Matplotlib figure, save it directly
    if isinstance(chord, plt.Figure):  
        chord.savefig(dest_plot, dpi=300, bbox_inches="tight")
    else:  
        # If `chord` is not a figure, get the current figure and save it
        plt.gcf().savefig(dest_plot, dpi=300, bbox_inches="tight")

    plt.close()

    print(f"Saved plot: {dest_plot}")


    


def start() -> None:
    import os

    # Load merged hs_names data
    dest = f"{checkpoint_dir}/adata_final_merged_raw_norm_hs_names.h5ad"
    if os.path.exists(dest):
        print("Load merged hs_names data...")
        print(dest)
        adata = sc.read_h5ad(dest)
    else:
        return
    
    ### Statistical analysis - Plotting
    # This method will retrieve interactions where the mean expression of the interacting partners
    # (proteins participating in the interaction) displays significant cell state specificity by
    # employing a random shuffling methodology.
    print(adata)
    if statistical_analysis:
        # Heatmaps
        # plot_heatmaps(adata)

        # # Lineages vs Other lineages interactions
        # plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Neu")
        # plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="MeV")
        # plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Imm")
        
        chord_diagram(adata=adata)
            
        
start()

print("\n********\n* DONE *\n********")
