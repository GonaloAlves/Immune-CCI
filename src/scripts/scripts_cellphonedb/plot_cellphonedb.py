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
    dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"
    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)


    # Generate heatmap
    clusterg = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        log1p_transform=log1p,
        figsize=(60, 60),
        linewidths=0.2,
        # annot = True,
        # annot_kws={"size": 20},
        col_cluster=True,  # prevent it from reordering automatically
        row_cluster=True,
        method='ward'
    )

    # Title customization
    suptitle = "Number of significant interactions (log1p transformed)" if log1p \
        else "Number of significant interactions"
    clusterg.figure.suptitle(suptitle, y=0.85, size=60)
    clusterg.ax_cbar.set_position((1.1, 0.0, .02, .3))  # Move color bar to the right
    clusterg.ax_cbar.tick_params(labelsize=50)  # Increase scale label font size

    # Adjust heatmap
    ax = clusterg.ax_heatmap
    ax.grid(False)

    # Increase font size of cluster labels
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=50)  # Adjust as needed
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=50, rotation=90)  # Rotate for readability
    
    # Save plot
    output_path = f"{cellphonedb_dir}/significant_interactions_final_merged_nona.png"
    print(f"Saving heatmap to: {output_path}")
    clusterg.savefig(output_path, bbox_inches="tight")
    plt.close()

def test_heatmap(adata: sc.AnnData, obs_key: str = None, category: str = None, remove_clusters: list = [], vmin: int = None, vmax: int = None):

    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap

    if obs_key is not None and category is None:
        raise ValueError("If obs_key is not None then category cannot be None!")
    elif obs_key is not None:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{category.replace('.', '_')}_nona.txt"
    else:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"

    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
    
    # Generate heatmap
    clusterg = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        return_tables=True,
        col_cluster=True, 
        row_cluster=True,
        method='ward'
    )
    
    # Extract the matrix
    count_matrix = clusterg["count_network"]
    
    custom_order = [
        "Imm.M0Like.1", "Imm.M0Like.2", "Imm.DAM.0", "Imm.MHCII.0", "Imm.Interferon.0", "Imm.PVM.0", "Imm.DAM.1", "Imm.Proliferative.0",
        "Neu.CSFcN.0", "Neu.Epend.0", "MeV.Epithelial.0", "MeV.Pericytes.0", "MeV.Endothelial.0", "MeV.SMC.0", "MeV.FibCollagen.1", "MeV.Fib.5", "MeV.Fib.4", 
        "MeV.FibProlif.0", "MeV.Endothelial.3", "MeV.FibCollagen.2", "MeV.VLMC.1", "MeV.FibLaminin.0", "MeV.VLMC.0", "MeV.Fib.3",
        "Imm.M0Like.0", "MeV.Endothelial.1", "MeV.Fib.2", "MeV.FibCollagen.0", "MeV.Fib.0", "MeV.Fib.1", "MeV.Endothelial.2", "MeV.FibCollagen.3"
    ]

    # Check if all your custom labels exist in the matrix
    missing = set(custom_order) - set(count_matrix.index)
    if missing:
        print("Warning: These cluster names are not in the matrix:", missing)

    # Remove clusters you want to exclude from the count_matrix
    if remove_clusters:
        print(f"Removing the following clusters: {remove_clusters}")
        count_matrix = count_matrix.drop(index=remove_clusters, columns=remove_clusters, errors='ignore')

    # Reorder the matrix rows and columns (if custom_order still applies)
    ordered_matrix = count_matrix.loc[custom_order, custom_order]
    print(f"Ordered Matrix:\n{ordered_matrix}")


    # print("######")
    # print(ordered_matrix.shape)
    # print("######")

    plt.figure(figsize=(60, 60))

    # Define your custom gradient colors: Blue → Light Beige → Purplish-Red
    custom_colors = ["#2166ac", "#ffead0", "#b2182b"]  # Replace with exact hex if needed

    # Create the colormap
    custom_cmap = LinearSegmentedColormap.from_list("custom_bluered", custom_colors, N=256)

    ax = sns.heatmap(
        ordered_matrix,
        annot=True,
        fmt=".0f",
        cmap=custom_cmap,
        linewidths=0.2,
        linecolor='gray',
        square=True,
        cbar_kws={"shrink": 0.8},
        xticklabels=True,
        yticklabels=True,
        annot_kws={"size": 30},
        vmin=vmin,
        vmax=vmax
    )

    # Title and axis customization
    plt.title(f"Number of Significant Interactions in {category}", fontsize=60, pad=60)
    ax.yaxis.set_ticks_position('right')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=40, rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=40, rotation=0)

    # Adjust colorbar ticksM0Like.1
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    cbar.ax.set_position([0.85, 0.2, 0.5, 0.3])  # [left, bottom, width, height]

    # Save the figure
    output_path = f"{cellphonedb_dir}/manual_heatmap_{category}.png"
    plt.savefig(output_path, bbox_inches="tight")
    plt.close()

    return ordered_matrix

def export_to_excel(df, output_path):
    """
    Exports the given DataFrame to an Excel file.
    
    Parameters:
    - df: pd.DataFrame
    - output_path: str, full path to save the Excel file
    """
    try:
        df.to_excel(output_path, index=False)
        print(f"DataFrame exported successfully to {output_path}")
    except Exception as e:
        print(f"Failed to export DataFrame: {e}")

def plot_heatmaps_fixed_order(adata: sc.AnnData, obs_key: str = None, category: str = None, log1p: bool = False) -> None:
    """
    Plots a heatmap using ktplotspy for a given AnnData object.

    Parameters
    ----------
    `adata` : sc.AnnData
        AnnData object used for plotting.
        
    `log1p` : bool, optional
        Whether to log-transform the number of significant interactions for better visualization.
    """


    
    if obs_key is not None and category is None:
        raise ValueError("If obs_key is not None then category cannot be None!")
    elif obs_key is not None:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{category.replace('.', '_')}.txt"

    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
    
    
    # Generate heatmap
    clusterg = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        # log1p_transform=log1p,
        # figsize=(60, 60),
        # linewidths=0.2,
        # annot = True,
        return_tables=True,
        col_cluster=True, 
        row_cluster=True,
        method='ward'
    )
    #print(clusterg)
    # Extract the matrix
    count_matrix = clusterg["count_network"]
    #print(count_matrix)

    custom_order = [
        "Imm.M0Like.1", "Imm.M0Like.2", "Imm.DAM.0", "Imm.MHCII.0", "Imm.Interferon.0", "MeV.ImmuneDoublets.0", "Imm.PVM.0", "Imm.DAM.1", "Imm.Proliferative.0",
        "Neu.CSFcN.0", "Neu.Epend.0", "MeV.Epithelial.0" ,"MeV.Pericytes.0", "MeV.Endothelial.0", "MeV.SMC.0", "MeV.FibCollagen.1", "MeV.Fib.5", "MeV.Fib.4", 
        "MeV.FibProlif.0", "MeV.Endothelial.3", "MeV.FibCollagen.2", "MeV.VLMC.1", "MeV.FibLaminin.0", "MeV.VLMC.0", "MeV.Fib.3", "MeV.FibUnknown.6", "MeV.LowQuality.0",
        "Imm.M0Like.0", "MeV.Endothelial.1", "MeV.Fib.2", "MeV.FibCollagen.0", "MeV.Fib.0", "MeV.Fib.1", "MeV.Endothelial.2", "MeV.FibCollagen.3"
    ]

    # Check if all your custom labels exist in the matrix
    missing = set(custom_order) - set(count_matrix.index)
    if missing:
        print("Warning: These cluster names are not in the matrix:", missing)

    # Reorder the matrix rows and columns
    ordered_matrix = count_matrix.loc[custom_order, custom_order]
    print(f"hey_{ordered_matrix}")

    # Generate heatmap
    clusterg2 = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        cell_types=ordered_matrix,  
        log1p_transform=log1p,
        figsize=(60, 60),
        linewidths=0.2,
        annot = True,
        col_cluster=False,  # prevent it from reordering automatically
        row_cluster=False
        # method='ward'
    )


    # Title customization
    suptitle = "Number of significant interactions (log1p transformed)" if log1p \
        else f"Number of significant interactions in {category}"
    clusterg2.figure.suptitle(suptitle, y=0.85, size=60)
    clusterg2.ax_cbar.set_position((1.1, 0.0, .02, .3))  # Move color bar to the right
    clusterg2.ax_cbar.tick_params(labelsize=50)  # Increase scale label font size

    # Adjust heatmap
    ax = clusterg2.ax_heatmap
    ax.grid(False)

    # Increase font size of cluster labels
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=40)  # Adjust as needed
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=40, rotation=90)  # Rotate for readability
    
    # Save plot
    output_path = f"{cellphonedb_dir}/significant_interactions_final_merged_nona_{category}.png"
    print(f"Saving heatmap to: {output_path}")
    clusterg2.savefig(output_path, bbox_inches="tight")
    plt.close()


def plot_heatmaps_auto_order(adata: sc.AnnData, obs_key: str = None, category: str = None, log1p: bool = False) -> None:
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
    dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"
    pvalues1 = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
    
    
    # Generate heatmap
    clusterg = kpy.plot_cpdb_heatmap(
        pvals=pvalues1,
        # log1p_transform=log1p,
        # figsize=(60, 60),
        # linewidths=0.2,
        # annot = True,
        return_tables=True,
        col_cluster=True, 
        row_cluster=True,
        method='ward'
    )
    print(clusterg)
    # Extract the matrix
    count_matrix = clusterg["count_network"]
    print(count_matrix)

    # Get the ordered list of clusters
    ordered_clusters = (
        count_matrix
        .index
        .tolist()
    )

    # Reorder the matrix rows and columns
    ordered_matrix = count_matrix.loc[ordered_clusters, ordered_clusters]
    print(ordered_matrix)

    if obs_key is not None and category is None:
        raise ValueError("If obs_key is not None then category cannot be None!")
    elif obs_key is not None:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona_{category.replace('.', '_')}.txt"
    print(f"loadin... {dest}")

    # Load p-values file (assuming a fixed filename)
    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)

    # Generate heatmap
    clusterg2 = kpy.plot_cpdb_heatmap(
        pvals=pvalues,
        cell_types=ordered_matrix,  
        log1p_transform=log1p,
        figsize=(60, 60),
        linewidths=0.2,
        annot = True,
        col_cluster=False,  # prevent it from reordering automatically
        row_cluster=False
        # method='ward'
    )


    # Title customization
    suptitle = "Number of significant interactions (log1p transformed)" if log1p \
        else "Number of significant interactions"
    clusterg2.figure.suptitle(suptitle, y=0.85, size=60)
    clusterg2.ax_cbar.set_position((1.1, 0.0, .02, .3))  # Move color bar to the right
    clusterg2.ax_cbar.tick_params(labelsize=50)  # Increase scale label font size

    # Adjust heatmap
    ax = clusterg2.ax_heatmap
    ax.grid(False)

    # Increase font size of cluster labels
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=40)  # Adjust as needed
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=40, rotation=90)  # Rotate for readability
    
    # Save plot
    output_path = f"{cellphonedb_dir}/significant_interactions_final_merged_nona_{category}.png"
    print(f"Saving heatmap to: {output_path}")
    clusterg2.savefig(output_path, bbox_inches="tight")
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
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{category.replace('.', '_')}_nona.txt"
    else:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"
        
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
        dest = cellphonedb_dir + f"/significant_interactions_final_merged_major_{category.replace('.', '_')}_nona.pdf"
    else:
        dest = cellphonedb_dir + "/significant_interactions_final_merged_major_nona.pdf"
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
    dest_pvalues = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"
    dest_means = f"{cellphonedb_dir}/statistical_analysis_means_final_merged_nona.txt"
    
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
    dest_plot = f"{cellphonedb_dir}/{lineage_prefix}_interactions_final_merged_nona.pdf"
    dotplot.save(dest_plot, dpi=300, limitsize=False, bbox_inches="tight")
    plt.close()

    print(f"Saved plot: {dest_plot}")


def chord_diagram(adata: sc.AnnData, lineage_prefix: str) -> None:

    # help(kpy.plot_cpdb_chord)
    
    cell_types = adata.obs["leiden_merge"].cat.categories.str.startswith(lineage_prefix)
    if not cell_types.any():
        raise ValueError(f"Lineage prefix '{lineage_prefix}' not found!")

    # Load interaction data (ignoring injury conditions)
    dest_pvalues = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"
    dest_means = f"{cellphonedb_dir}/statistical_analysis_means_final_merged_nona.txt"

    pvalues = pd.read_csv(dest_pvalues, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    print(pvalues)
    means = pd.read_csv(dest_means, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    print(means)

    all_cell_types = adata.obs["leiden_merge"].cat.categories.to_list()
    # Use all cell types for both groups to compare everything against everything
    cell_types1 = [cell for cell in all_cell_types if cell.startswith(lineage_prefix)]
    cell_types2 = [cell for cell in all_cell_types if not cell.startswith(lineage_prefix)]

    cell_types1 = "|".join(cell_types1)
    cell_types2 = "|".join(cell_types2)

    deconvoluted_file = f"{cellphonedb_dir}/statistical_analysis_deconvoluted_final_merged_nona.txt"
    deconvoluted_data = pd.read_csv(deconvoluted_file, sep="\t")

##################
    #Debugging

    # # METHOD1 
    # # Check cell types in adata
    # adata_cell_types = set(adata.obs["leiden_merge"].unique())

    # valid_cell_types = adata_cell_types.intersection(deconvoluted_data)

    # adata = adata[adata.obs["leiden_merge"].isin(valid_cell_types)]
    # deconvoluted_data = deconvoluted_data[list(valid_cell_types)]
    
    #####
    
    # # Find mismatches
    # only_in_adata = adata_cell_types - deconvoluted_cell_types
    # only_in_deconvoluted = deconvoluted_cell_types - adata_cell_types

    # print(f"Cell types only in adata: {only_in_adata}")
    # print(f"Cell types only in adata: {adata_cell_types}")
    # print(f"Cell types only in deconvoluted data: {only_in_deconvoluted}")
    # print(f"Cell types only in adata: {deconvoluted_cell_types}")


    # # METHOD2 Ensure all cell type names are lowercase and without spaces/dots
    # adata.obs["leiden_merge"] = adata.obs["leiden_merge"].str.replace(" ", "_").str.lower()
    # deconvoluted_data["interacting_pair"] = deconvoluted_data["interacting_pair"].str.replace(" ", "_").str.lower()

    # # METHOD3 Keep only matching cell types
    # valid_cell_types = adata_cell_types.intersection(deconvoluted_cell_types)
    # print(valid_cell_types)

    # print(len(complex_id), tmpdf.shape[0])

    # print("Means shape:", means.shape)
    # print("P-values shape:", pvalues.shape)
    # print("Deconvoluted data shape:", deconvoluted_data.shape)

    # adata = adata[adata.obs["leiden_merge"].isin(valid_cell_types)]
    # deconvoluted_data = deconvoluted_data[deconvoluted_data["interacting_pair"].isin(valid_cell_types)]
# ##################



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
        # link_kwargs={"direction": 1, "allow_twist": True, "r1": 95, "r2": 90},
        # sector_text_kwargs={"color": "black", "size": 12, "r": 105, "adjust_rotation": True},
        # legend_kwargs={"loc": "center", "bbox_to_anchor": (1, 1), "fontsize": 8},
        max_size=3,
        highlight_size=1,
        standard_scale=True
    )

    
    # Save plot
    dest_plot = f"{cellphonedb_dir}/chord_all_interactions_final_merged_nona.pdf"

    # If `chord` is a Matplotlib figure, save it directly
    if isinstance(chord, plt.Figure):  
        chord.savefig(dest_plot, dpi=300, bbox_inches="tight")
    else:  
        # If `chord` is not a figure, get the current figure and save it
        plt.gcf().savefig(dest_plot, dpi=300, bbox_inches="tight")

    plt.close()

    print(f"Saved plot: {dest_plot}")


def schord_diagram(adata: sc.AnnData) -> None:

    # help(kpy.plot_cpdb_chord)
    
    
    # Load interaction data (ignoring injury conditions)
    dest_pvalues = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_nona.txt"
    dest_means = f"{cellphonedb_dir}/statistical_analysis_means_final_merged_nona.txt"

    pvalues = pd.read_csv(dest_pvalues, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    print(pvalues)
    means = pd.read_csv(dest_means, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    print(means)

    all_cell_types = adata.obs["leiden_merge"].cat.categories.to_list()
    print(adata.obs["leiden_merge"].cat.categories.to_list())
    
    # Use all cell types for both groups to compare everything against everything
    
    cell_types1 = "|".join(all_cell_types)
    cell_types2 = "|".join(all_cell_types)


    deconvoluted_file = f"{cellphonedb_dir}/statistical_analysis_deconvoluted_final_merged_nona.txt"
    deconvoluted_data = pd.read_csv(deconvoluted_file, sep="\t")


    # METHOD1 
    # # Check cell types in adata
    # adata_cell_types = set(adata.obs["leiden_merge"].unique())

    # valid_cell_types = adata_cell_types.intersection(deconvoluted_data)

    # adata = adata[adata.obs["leiden_merge"].isin(valid_cell_types)]
    # deconvoluted_data = deconvoluted_data[list(valid_cell_types)]
    


    chord = kpy.plot_cpdb_chord(
        adata=adata,
        cell_type1=".",
        cell_type2=".",
        means=means,
        pvals=pvalues,
        deconvoluted= deconvoluted_data,
        celltype_key="leiden_merge",
        figsize=(40, 30),
        title="Interactions between all cells",
        # link_kwargs={"direction": 1, "allow_twist": True, "r1": 95, "r2": 90},
        # sector_text_kwargs={"color": "black", "size": 12, "r": 105, "adjust_rotation": True},
        # legend_kwargs={"loc": "center", "bbox_to_anchor": (1, 1), "fontsize": 8},
        max_size=3,
        highlight_size=1,
        standard_scale=True
    )

    
    # Save plot
    dest_plot = f"{cellphonedb_dir}/chord_all_interactions_final_merged_nona.pdf"

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
    dest = f"{checkpoint_dir}/adata_final_merged_raw_norm_hs_names_nona.h5ad"
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
        
        # plot_heatmaps_fixed_order(adata, obs_key="injury_day", category="sham_15")
        # plot_heatmaps_fixed_order(adata, obs_key="injury_day", category="injured_15")
        # plot_heatmaps_fixed_order(adata, obs_key="injury_day", category="injured_60")
        
        remove_clusters = ["MeV.ImmuneDoublets.0", "MeV.FibUnknown.6", "MeV.LowQuality.0"]
        test_heatmap(adata)
        test = test_heatmap(adata, obs_key="injury_day", category="injured_15", remove_clusters=remove_clusters, vmin = 0, vmax = 95)
        test_heatmap(adata, obs_key="injury_day", category="injured_60", remove_clusters=remove_clusters, vmin = 0, vmax = 95)
        test_heatmap(adata, obs_key="injury_day", category="uninjured", remove_clusters=remove_clusters, vmin = 0, vmax = 95)

        export_to_excel(test, "test_simplified.xlsx")

        # # Lineages vs Other lineages interactions
        # plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Neu")
        # plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="MeV")
        # plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Imm")
        
        #chord_diagram(adata=adata, lineage_prefix="Imm")
        # schord_diagram(adata=adata)
            
        
start()

print("\n********\n* DONE *\n********")
