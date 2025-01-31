import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp


proportion_colors = {
    'all': '#899499',
    'uninjured_0': '#428bca',
    'sham_15': '#5cb85c',
    'injured_15': '#ff9922',
    'injured_60': '#d9534f',
    'uninjured': '#428bca',
    'injured': '#ff9922',
    'sham': '#5cb85c',
    'no_sci': '#428bca',
    'sci': '#ff9922'
}

def load_data(file_path):
    """
    Load an AnnData object from an .h5ad file.
    """
    print(f"Loading data from: {file_path}")
    return sc.read_h5ad(file_path)

def get_sample_groups(dataset_type, condition):
    """
    Get the list of sample groups based on the condition type.
    """
    sample_groups = {
        'injury_day': ['uninjured_0', 'sham_15', 'injured_15', 'injured_60'],
        'injury': ['uninjured', 'sham', 'injured'],
        'injury_grouped': ['no_sci', 'sci']
    }
    
    samples = sample_groups.get(condition, [])
    if 'uinj' in dataset_type:
        samples = [s for s in samples if not re.match(r'^injured', s, re.IGNORECASE)]
    
    return samples

def group_injury_conditions(adata, original_column='injury', new_column='injury_grouped'):
    """
    Create a new column that groups 'uninjured' and 'sham' into 'uninjured_combined'.

    Parameters:
    - adata (AnnData): The AnnData object with metadata.
    - original_column (str): The column to map from (e.g., 'injury').
    - new_column (str): The new column to store grouped conditions.

    Returns:
    - AnnData with the new column added.
    """
    mapping = {
        'uninjured': 'no_sci',
        'sham': 'no_sci',
        'injured': 'sci'  
    }

    # Create the new column in adata.obs
    adata.obs[new_column] = adata.obs[original_column].replace(mapping)

    print(f"New column '{new_column}' created in adata.obs.")

    print(f"{new_column} data: {adata.obs[new_column]}")

    return adata


def calculate_cell_fractions(adata, dataset_type, key, condition):
    """
    Compute cell fractions per cluster for a given condition.
    """
    samples = get_sample_groups(dataset_type, condition)

    print(samples)
    print("Keys:", list(adata.obs.columns))
    
    adata_filtered = adata[adata.obs[condition].isin(samples), :].copy()
    cluster_categories = adata_filtered.obs[key].cat.categories
    sample_counts = adata_filtered.obs[condition].value_counts().loc[samples].copy()
    
    # Initialize table to store fractions
    table_frac = np.zeros((len(cluster_categories) + 2, len(samples) + 1))
    
    total_cluster = adata_filtered.obs[key].value_counts().loc[cluster_categories].values.tolist()
    total_sample = sample_counts.values.tolist()
    
    for i, cluster in enumerate(cluster_categories):
        for j, sample in enumerate(samples):
            n_cells = adata_filtered[(adata_filtered.obs[key] == cluster) &
                                     (adata_filtered.obs[condition] == sample)].n_obs
            table_frac[i, j] = n_cells / total_cluster[i]
    
    # Calculate sample fractions in total
    row_total = len(total_cluster)
    for j in range(len(samples)):
        table_frac[row_total, j] = adata_filtered[adata_filtered.obs[condition] == samples[j]].n_obs / adata_filtered.n_obs
        table_frac[row_total + 1, j] = total_sample[j]
    
    table_frac[row_total, -1] = adata_filtered.n_obs / adata_filtered.n_obs
    table_frac[row_total + 1, -1] = adata_filtered.n_obs
    
    row_names = list(cluster_categories) + ['frac sample in total', 'total in sample']
    col_names = samples + ['total in cluster']
    
    return pd.DataFrame(table_frac, index=row_names, columns=col_names)

def save_fractions_table(table, output_dir, dataset, clusters, condition):
    """
    Save the computed cell fractions to a text file.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{dataset}_sample_frac_{clusters}_{condition}.txt")
    table.to_csv(output_file, sep='\t')
    print(f"Saved fractions table to: {output_file}")

def plot_stacked_bar(data, dataset_name, output_dir, alternate_names=None):
    """
    Generate and save a stacked bar plot for cell fractions.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{dataset_name}_stacked_plot.pdf")
    
    clusters = range(len(data.index))
    bottom = np.zeros(len(clusters))
    colors = [proportion_colors[c] for c in data.columns]
    
    with PdfPages(output_file, keep_empty=True) as pdf:
        fig, ax = plt.subplots(figsize=(len(data.index) * 0.5, 7))
        ax.grid(False)
        
        for i, col in enumerate(data.columns):
            ax.bar(clusters, height=data.iloc[:, i].values, bottom=bottom, color= colors[i], width=0.6, label=col)
            bottom += data.iloc[:, i].values
        
        ax.set_xticks(clusters)
        ax.set_xticklabels(alternate_names or data.index, rotation=45, ha='right', rotation_mode='anchor')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        pdf.savefig(dpi=300, bbox_inches='tight')
        plt.close(fig)
    print(f"Saved plot: {output_file}")

def start_analysis(input_file, output_dir):
    """
    Load data, compute cell fractions, and generate stacked bar plots.

    Parameters:
    - input_file (str): Path to the .h5ad file.
    - output_dir (str): Directory to save output files.

    Returns:
    - None
    """

    # Load AnnData object
    adata = load_data(input_file)

    # Create the new grouped injury column
    adata = group_injury_conditions(adata, original_column='injury', new_column='injury_grouped')

    clusters_key = 'leiden_fusion'

    # Process each condition type
    for condition in ['injury_day', 'injury', 'injury_grouped']:  
        print(f"Processing condition: {condition}...")

        table_frac = calculate_cell_fractions(adata, 'Immune', clusters_key, condition)
        
        # Save fraction table
        output_file = os.path.join(output_dir, f"Immune_sample_frac_{clusters_key}_{condition}.txt")
        print(f"Saving table: {output_file}")
        table_frac.to_csv(output_file, sep='\t')

        # Store in adata.uns for later analysis
        adata.uns[f'sample_frac_{clusters_key}_{condition}'] = table_frac
        
        # Get expected frequencies
        expected_freq = table_frac.iloc[-2, :-1].values

        # Plot stacked barplots
        plot_data = table_frac.iloc[:-2, :-1].copy()
        plot_data.loc['Expected', :] = expected_freq
        alternate_names = list(table_frac.index[:-2]) + ['Expected']

        plot_stacked_bar(plot_data, 
                         dataset_name=f'Immune_final_{clusters_key}_{condition}', 
                         output_dir=output_dir, 
                         alternate_names=alternate_names)    
            
    # adata.write_h5ad(input_file, compression='gzip')
    # print("Analysis completed.")

if __name__ == "__main__":
    try:
        mp.set_start_method('spawn', force=True)
        input_file = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy_copy.h5ad"
        output_dir = "/home/makowlg/Documents/Immune-CCI/src/fractions_related/Immune"
        start_analysis(input_file, output_dir)
        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
