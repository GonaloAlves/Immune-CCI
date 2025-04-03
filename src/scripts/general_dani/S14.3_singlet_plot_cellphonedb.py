# Perform CellPhoneDB analysis - LeonorSaude10x
#
# Follwing these tutorials
# https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T01_Method1.ipynb
# https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T01_Method2.ipynb
#
# Daniel Ribeiro, 2023
import gc

import pandas as pd
import scanpy as sc

statistical_analysis = True
deg_analysis = True


def plot_heatmaps(adata: sc.AnnData,
                  obs_key: str = None,
                  category: str = None,
                  log1p: bool = False
                  ) -> None:
    """
    Parameters
    ----------
    `adata` : sc.AnnData
        AnnData object. Passed to ktplotspy plot function.
        
    `obs_key` : str
        Key to retreive available categories. Not fully implemented, except to evaluate whether \
            somenone asks for a category without providing a obs_key.
    
    `category` : str
        Category name to be found in .txt file name.
        
    `log1p` : bool
        Wether to log the number of significant categories. May improve heatmap vizualization. \
        Passed to ktplotspy.
    """
    import matplotlib
    matplotlib.use('AGG')   # Cairo doesnt work here because draw_gouraud_triangle is not implemented for png
    import matplotlib.pyplot as plt
    import ktplotspy as kpy
    from src.globals import cellphonedb_dir
    
    if obs_key is not None and category is None:
        raise ValueError("If obs_key is not None then category cannot be None!")
    elif obs_key is not None:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{category.replace('.', '_')}.txt"
    else:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged.txt"
    
    # Plot heatmap
    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
    clusterg = kpy.plot_cpdb_heatmap(
        adata=adata,
        pvals=pvalues,
        log1p_transform=log1p,
        celltype_key="leiden_merge",
        figsize=(60, 60),
        linewidths=0.2)
    suptitle = "Number of significant interactions (log1p transformed)" if log1p \
        else "Number of significant interactions"
    clusterg.figure.suptitle(suptitle, y=0.85, size=50)
    clusterg.ax_cbar.set_position((1.02, 0.0, .02, .3))
    ax = clusterg.ax_heatmap
    ax.grid(False)
    plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False)
    
    print(dest)
    clusterg.savefig(dest, bbox_inches="tight")
    plt.close()


def plot_heatmaps_major_cells(obs_key: str = None,
                              category: str | None = None,
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
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from src.globals import cellphonedb_dir
    from src.helper_functions import get_cell_types
    
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


def plot_cluster_72_interactions(adata: sc.AnnData,
                                 genes: list = None,
                                 suffix: str = None
                                 ) -> None:
    import matplotlib
    matplotlib.use('AGG')   # Cairo doesnt work here because draw_gouraud_triangle is not implemented for png
    import matplotlib.pyplot as plt
    import ktplotspy as kpy
    from src.globals import cellphonedb_dir
    
    # Cluster 72 data only makes sense for inured_15 and injured_60
    # Cluster 72 interaction - injured_15
    # Get cell types
    adata2 = adata[adata.obs["injury_day"] == "injured.15"].copy()
    dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_injured_15.txt"
    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    dest = f"{cellphonedb_dir}/statistical_analysis_means_final_merged_injured_15.txt"
    means = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    cell_types2 = adata2.obs["leiden_merge"].cat.categories.to_list()
    cell_types2 = "|".join(cell_types2)
    cell_type1 = "Neu.72"
    dotplot = kpy.plot_cpdb(
        adata=adata2,
        cell_type1="Neu.72",
        cell_type2=cell_type1,
        means=means,
        pvals=pvalues,
        celltype_key="leiden_merge",
        #default_style=False,
        genes=genes,
        figsize=(50, 25),
        title=f"Interactions between {cell_type1} and rest",
        max_size=6,
        highlight_size=2,
        standard_scale=True)
    if suffix is None:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_final_merged_injured_15.pdf"
    else:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_final_merged_injured_15_{suffix}.pdf"
    dotplot.save(dest, dpi=300, limitsize=False, bbox_inches="tight")
    plt.close()
    gc.collect()
    # Cluster 72 interaction - injured_15
    # Get cell types except other neurons
    cell_types2 = adata2.obs["leiden_merge"].cat.categories.to_list()
    cell_types2 = [cell for cell in cell_types2 if not cell.startswith("Neu")]
    cell_types2 = "|".join(cell_types2)
    dotplot = kpy.plot_cpdb(
        adata=adata2,
        cell_type1=cell_type1,
        cell_type2=cell_types2,
        means=means,
        pvals=pvalues,
        celltype_key="leiden_merge",
        #default_style=False,
        genes=genes,
        figsize=(50, 25),
        title=f"Interactions between {cell_type1} and rest",
        max_size=6,
        highlight_size=2,
        standard_scale=True)
    if suffix is None:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_noNeu_final_merged_injured_15.pdf"
    else:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_noNeu_final_merged_injured_15_{suffix}.pdf"
    dotplot.save(dest, dpi=300, limitsize=False, bbox_inches="tight")
    plt.close()
    gc.collect()
    
    # Cluster 72 interaction - injured_60
    # Get cell types
    adata2 = adata[adata.obs["injury_day"] == "injured.60"].copy()
    dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_injured_60.txt"
    pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    dest = f"{cellphonedb_dir}/statistical_analysis_means_final_merged_injured_60.txt"
    means = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    cell_types2 = adata2.obs["leiden_merge"].cat.categories.to_list()
    cell_types2 = "|".join(cell_types2)
    dotplot = kpy.plot_cpdb(
        adata=adata2,
        cell_type1=cell_type1,
        cell_type2=cell_types2,
        means=means,
        pvals=pvalues,
        celltype_key="leiden_merge",
        #default_style=False,
        genes=genes,
        figsize=(50, 25),
        title=f"Interactions between {cell_type1} and rest",
        max_size=6,
        highlight_size=2,
        standard_scale=True)
    if suffix is None:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_final_merged_injured_60.pdf"
    else:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_final_merged_injured_60_{suffix}.pdf"
    dotplot.save(dest, dpi=300, limitsize=False, bbox_inches="tight")
    plt.close()
    gc.collect()
    # Cluster 72 interaction - injured_60
    # Get cell types except other neurons
    cell_types2 = adata2.obs["leiden_merge"].cat.categories.to_list()
    cell_types2 = [cell for cell in cell_types2 if not cell.startswith("Neu")]
    cell_types2 = "|".join(cell_types2)
    dotplot = kpy.plot_cpdb(
        adata=adata2,
        cell_type1=cell_type1,
        cell_type2=cell_types2,
        means=means,
        pvals=pvalues,
        celltype_key="leiden_merge",
        #default_style=False,
        genes=genes,
        figsize=(50, 25),
        title=f"Interactions between {cell_type1} and rest",
        max_size=6,
        highlight_size=2,
        standard_scale=True)
    if suffix is None:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_noNeu_final_merged_injured_60.pdf"
    else:
        dest = cellphonedb_dir + f"/{cell_type1}_interactions_noNeu_final_merged_injured_60_{suffix}.pdf"
    dotplot.save(dest, dpi=300, limitsize=False, bbox_inches="tight")
    plt.close()
    gc.collect()


def plot_lineage_vs_other_interactions(adata: sc.AnnData,
                                       lineage_prefix: str) -> None:
    import matplotlib
    matplotlib.use('AGG')   # Cairo doesnt work here because draw_gouraud_triangle is not implemented for png
    import matplotlib.pyplot as plt
    import ktplotspy as kpy
    from src.globals import cellphonedb_dir
    
    injury_day = ["sham.15", "injured.15", "injured.60"]
    
    # Check if lineage exits
    cell_types = adata.obs["leiden_merge"].cat.categories.str.startswith(lineage_prefix)
    if not cell_types.any():
        raise ValueError(f"Lineage prefix {lineage_prefix} not found!")
    
    for inj in injury_day:
        adata2 = adata[adata.obs["injury_day"] == inj].copy()
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{inj.replace('.', '_')}.txt"
        pvalues = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
        dest = f"{cellphonedb_dir}/statistical_analysis_means_final_merged_{inj.replace('.', '_')}.txt"
        means = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
        # Select Lineage
        cell_types1 = adata2.obs["leiden_merge"].cat.categories.to_list()
        cell_types1 = [cell for cell in cell_types1 if cell.startswith(lineage_prefix)]
        cell_types1 = "|".join(cell_types1)
        # Select Other lineage
        cell_types2 = adata2.obs["leiden_merge"].cat.categories.to_list()
        cell_types2 = [cell for cell in cell_types2 if not cell.startswith(lineage_prefix)]
        cell_types2 = "|".join(cell_types2)
        # Plot
        dotplot = kpy.plot_cpdb(
            adata=adata2,
            cell_type1=cell_types1,
            cell_type2=cell_types2,
            means=means,
            pvals=pvalues,
            celltype_key="leiden_merge",
            figsize=(60, 50),
            title=f"Interactions between '{lineage_prefix}' and rest",
            max_size=3,
            highlight_size=1,
            standard_scale=True)
        dest = cellphonedb_dir + f"/{lineage_prefix}_interactions_final_merged_{inj.replace('.', '_')}.pdf"
        dotplot.save(dest, dpi=300, limitsize=False, bbox_inches="tight")
        plt.close()


def start() -> None:
    import os
    import src.globals   # noqa:F401
    from src.globals import checkpoint_dir

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
    if statistical_analysis:
        # Heatmaps
        plot_heatmaps(adata=adata)
        plot_heatmaps(adata=adata, obs_key="injury_day", category="sham_15")
        plot_heatmaps(adata=adata, obs_key="injury_day", category="injured_15")
        plot_heatmaps(adata=adata, obs_key="injury_day", category="injured_60")
        
        plot_heatmaps_major_cells(obs_key="injury_day", category="sham_15", log1p=True)
        plot_heatmaps_major_cells(obs_key="injury_day", category="injured_15", log1p=True)
        plot_heatmaps_major_cells(obs_key="injury_day", category="injured_60", log1p=True)
        
        # Calculate number of
        
        # Cluster 73 interaction
        #plot_cluster_73_interactions(adata=adata)
        #plot_cluster_73_interactions(adata=adata,
        #                             genes=["IL6", "IL6R", "IL6ST", "LIF", "LIFR", "CFTR", "TNF"],
        #                             suffix="genes")
        # Lineages vs Other lineages interactions
        #plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Neu")
        #plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Oli")
        #plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Ast")
        #plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="MeV")
        #plot_lineage_vs_other_interactions(adata=adata, lineage_prefix="Imm")
            
        
start()

print("\n********\n* DONE *\n********")
