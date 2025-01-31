# Calculate cell fractions in clusters - LeonorSaude10x
# 1. Calculate if some cluster is enriched for a particular sample
# 2. Calculate if some cluster is enriched for a particluar cell type seen by scMCA.
# 3. Plot fractions per cluster
#
#
# Daniel Ribeiro, 2023, ft Gonçalo Alves, 2025
import os
import re
import numpy as np
import pandas as pd

import scanpy as sc
import matplotlib
matplotlib.use('cairo')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

checkpoint_dir = "/home/makowlg/Documents/Immune-CCI/h5ad_files"
fractions_dir = "/home/makowlg/Documents/Immune-CCI/src/fractions_related/Immune"
data: dict[str, sc.AnnData] = {}

proportion_colors = {
    'all': '#899499',
    'uninjured_0': '#428bca',
    'sham_15': '#5cb85c',
    'injured_15': '#ff9922',
    'injured_60': '#d9534f',
    'uninjured_0_central': '#428bca',
    'sham_15_rostral': '#5cb85c',
    'sham_15_caudal': '#addbad',
    'injured_15_rostral': '#ff9922',
    'injured_15_caudal': '#f7bf81',
    'injured_60_rostral': '#d9534f',
    'injured_60_caudal': '#eca9a7',
    'caudal': '#b3bfa1',
    'rostral': '#d2a295',
    'central': '#428bca',
    'injured_caudal': '#9b636b',
    'injured_rostral': '#d2a295',
    'sham_caudal': '#5cb85c',
    'sham_rostral': '#addbad',
    'uninjured_central': '#428bca',
}

injury_day = ['uninjured_0', 'sham_15', 'injured_15', 'injured_60']
injury_condition = ['uninjured_0_central', 'sham_15_rostral', 'sham_15_caudal', 'injured_15_rostral', 'injured_15_caudal', 'injured_60_rostral', 'injured_60_caudal']
injury_region = ['uninjured_central', 'sham_rostral', 'sham_caudal', 'injured_rostral', 'injured_caudal']
collection_region = ['central', 'rostral', 'caudal']
injury_region_no_central = ['sham_rostral', 'sham_caudal', 'injured_rostral', 'injured_caudal']
collection_region_no_central = ['rostral', 'caudal']


def sample_factions(adata: sc.AnnData,
                    dataset_type: str,
                    key: str,
                    condition: str = 'injury_condition') -> pd.DataFrame:
    """
    adata: `sc.AnnData`
        Dataset containing cluster information
    key: `str`
        .obs key to look at
    condition: `str`
        String to choose from injury_day or injury_condition list of samples
    """


    samples = []
    prog = re.compile(r'^injured', flags=re.IGNORECASE)
    if condition == 'injury_condition':
        samples = list(injury_condition)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_day':
        samples = list(injury_day)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_region':
        samples = list(injury_region)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'collection_region':
        samples = list(collection_region)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'injury_region_no_central':
        condition = 'injury_region'
        samples = list(injury_region_no_central)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    elif condition == 'collection_region_no_central':
        condition = 'collection_region'
        samples = list(collection_region_no_central)
        if dataset_type.find('uinj') != -1:
            samples = [s for s in samples if prog.search(s) is None]
    else:
        raise ValueError("'condition' must be either a string 'injury_condition', 'injury_day', 'injury_region', 'collection_region'.")
    
    # Subset adata
    print("###########################Keys:", list(adata.obs.columns))
    adata2 = adata[adata.obs[condition].isin(samples), :].copy()
    # Initiate matrix
    n_cols = len(samples) + 1
    n_rows = len(adata2.obs[key].cat.categories) + 2
    table_frac = np.zeros((n_rows, n_cols))
    total_cluster = adata2.obs[key].value_counts().loc[adata2.obs[key].cat.categories].values.tolist()
    # Reorder samples logically
    sample_counts = adata2.obs[condition].value_counts().loc[samples].copy()
    total_sample = sample_counts.values.tolist()
    # Calculate fractions
    for i in range(len(total_cluster)):
        for j in range(len(sample_counts.index)):
            # Cells of cluster i that belong to sample j
            n_cells = adata2[(adata2.obs[key] == adata2.obs[key].cat.categories[i]) & \
                             (adata2.obs[condition] == sample_counts.index[j])].n_obs
            table_frac[i, j] = n_cells / total_cluster[i]
    # Calculate sample fractions in total and fill totals
    row_i = len(total_cluster)    # fraction in total
    for j in range(len(sample_counts.index)):
        table_frac[row_i, j] = adata2[adata2.obs[condition] == sample_counts.index[j]].n_obs / adata2.n_obs
        table_frac[row_i + 1, j] = total_sample[j]
    # Fill cluster totals
    col_j = len(total_sample)    # totals per cluster
    for i in range(len(total_cluster)):
        table_frac[i, col_j] = total_cluster[i]
    # Total cells in dataset
    table_frac[n_rows - 1, n_cols - 1] = adata2.n_obs
    table_frac[n_rows - 2, n_cols - 1] = adata2.n_obs / adata2.n_obs

    # Export table
    row_names = adata2.obs[key].cat.categories.to_list()
    row_names.append('frac sample in total')
    row_names.append('total in sample')
    col_names = samples
    col_names.append('total in cluster')
    table_frac = pd.DataFrame(table_frac, index=row_names, columns=col_names)

    return table_frac


def plot_stacked_plots(data: pd.DataFrame,
                       dataset_name: str,
                       alternate_sample_names=None
                       ) -> None:
    
    
    clusters = [i for i in range(len(data.index))]
    colors = [proportion_colors[c] for c in data.columns]
    bottom = np.zeros(len(clusters))

    with PdfPages(fractions_dir + f"/{dataset_name}_stacked_plot.pdf", keep_empty=True) as pdf:
        fig, ax = plt.subplots(figsize=(len(data.index) * 0.5, 7))
        ax.grid(False)
        for i, b in enumerate(data.columns):
            ax.bar(clusters,
                   height=data.iloc[:, i].to_numpy().flatten(),
                   bottom=bottom,
                   color=colors[i],
                   width=0.6,
                   yerr=0,
                   label=b)
            # Next iteration bars will be plotted on top of the previous one
            bottom += data.iloc[:, i].to_numpy().flatten()

        ax.set_xticks(clusters)
        if alternate_sample_names is not None:
            ax.set_xticklabels(alternate_sample_names, rotation=45, ha='right', rotation_mode='anchor')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        pdf.savefig(dpi=300, bbox_inches='tight')
        plt.close(fig)


def start(n_proc=None) -> None:

    # Load ranked data
    #for d in list(datasets_divided.keys()):
    #    datasets_divided[f'{d}_uinj'] = datasets_divided[d]
   
    d = 'Immune'
    dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked_copy_copy.h5ad"
    if os.path.exists(dest):
        print("Load gene rank data...")
        print(dest)
        data[d] = sc.read_h5ad(dest)

    
    # Calculate fractions and totals
    
    clusters = f'leiden_fusion'
    for condition in ['injury_condition', 'injury_day', 'injury_region', 'collection_region',
                        'injury_region_no_central', 'collection_region_no_central']:
        print(f"Calculate statistics for {condition}...")
        # Save fraction
        table_frac = sample_factions(data[d], dataset_type=d, key=clusters, condition=condition)
        dest = f"{fractions_dir}/{d}_final_sample_frac_{clusters}_{condition}.txt"
        print(f"Output: {dest}")
        table_frac.to_csv(dest, sep='\t')
        # Update AnnData
        data[d].uns[f'sample_frac_{clusters}_{condition}'] = table_frac
        
        expected_freq = table_frac.iloc[-2, :-1].values  # Expected freqeuncies observed in data set

        
        # Plot stacked barplots
        print("\nPlot stacked barplots...")
        plot_data = data[d].uns[f'sample_frac_{clusters}_{condition}'].iloc[:-2, :-1].copy()
        plot_data.loc['Expected', :] = expected_freq
        alternate_names = data[d].uns[f'sample_frac_{clusters}_{condition}'].index[:-2].to_list()
        alternate_names.append('Expected')
        plot_stacked_plots(plot_data,
                            dataset_name=f'{d}_final_{clusters}_{condition}',
                            alternate_sample_names=alternate_names)
        # With dendrogram order
        plot_data = data[d].uns[f'sample_frac_{clusters}_{condition}'].iloc[:-2, :-1].copy()
        plot_data = plot_data.loc[data[d].uns[f'dendrogram_{clusters}']['dendrogram_info']['ivl'].tolist(), :]
        plot_data.loc['Expected', :] = expected_freq
        alternate_names = data[d].uns[f'dendrogram_{clusters}']['dendrogram_info']['ivl'].tolist()
        alternate_names.append('Expected')
        plot_stacked_plots(plot_data,
                            dataset_name=f'{d}_final_{clusters}_{condition}_dendrogram',
                            alternate_sample_names=alternate_names)
    
    # Save AnnData
    print("Save AnnData...")
    dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked_copy_copy.h5ad"
    data[d].write_h5ad(dest, compression='gzip')


# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=8)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
