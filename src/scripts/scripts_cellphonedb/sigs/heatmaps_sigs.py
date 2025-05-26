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
cellphonedb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/filtered_pvalues"
cellphonedb_dir_plots = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/plots/heatmaps/sigs"

def test_heatmap(obs_key: str = None, category: str = None, vmin: int = None, vmax: int = None, restrictions = None):

    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from matplotlib.colors import LinearSegmentedColormap

    if obs_key is not None and category is None:
        raise ValueError("If obs_key is not None then category cannot be None!")
    elif obs_key is not None:
        dest = f"{cellphonedb_dir}/statistical_analysis_pvalues_final_merged_{category.replace('.', '_')}_nona_{restrictions}.txt"
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
        "Imm.M0Like.1", "Imm.DAM.0", "Imm.Interferon.0", "Imm.PVM.0", "Imm.DAM.1",
        "Neu.CSFcN.0", "Neu.Epend.0", "MeV.Pericytes.0", "MeV.Endothelial.0", "MeV.FibCollagen.1", "MeV.Fib.5", "MeV.Fib.4", 
        "MeV.FibCollagen.2", "MeV.Endothelial.1", "MeV.Endothelial.2", "MeV.FibCollagen.3"
    ]

    relevant_clusters = ["Imm.DAM.0", "Imm.Interferon.0","Imm.PVM.0", "Imm.DAM.1","Neu.Epend.0", "MeV.Pericytes.0", "MeV.Endothelial.0","MeV.Endothelial.1", "MeV.Endothelial.2"]


    # Check if all your custom labels exist in the matrix
    missing = set(relevant_clusters) - set(count_matrix.index)
    if missing:
        print("Warning: These cluster names are not in the matrix:", missing)

    # Reorder the matrix rows and columns (if custom_order still applies)
    ordered_matrix = count_matrix.loc[relevant_clusters, relevant_clusters]
    print(f"Ordered Matrix:\n{ordered_matrix}")

    #show only one part
    # mask = np.triu(np.ones(ordered_matrix.shape, dtype=bool), k=1)
    mask = np.tril(np.ones(ordered_matrix.shape, dtype=bool), k=-1) 
    
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
        mask=mask,
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

    # Move x-tick labels to top
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=40, rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=40, rotation=0)
    ax.yaxis.set_ticks_position('right')

    # Add title at the bottom manually
    plt.figtext(0.5, 0.01, f"Number of Significant Interactions in {category}", ha='center', fontsize=60)

    # Adjust colorbar ticksM0Like.1
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=40)
    cbar.ax.set_position([0.85, 0.2, 0.5, 0.3])  # [left, bottom, width, height]

    # Save the figure
    output_path = f"{cellphonedb_dir_plots}/manual_heatmap_{category}_{restrictions}.png"
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
        
        # remove_clusters = ["MeV.ImmuneDoublets.0", "MeV.FibUnknown.6", "MeV.LowQuality.0"]

        test_heatmap(obs_key="injury_day", category="injured_15", vmin = 0, vmax = 90, restrictions="two")
        test_heatmap(obs_key="injury_day", category="injured_60", vmin = 0, vmax = 90, restrictions="two")
        test_heatmap(obs_key="injury_day", category="uninjured", vmin = 0, vmax = 90, restrictions="two")

        # export_to_excel(test, "test_simplified.xlsx")

            
        
start()

print("\n********\n* DONE *\n********")