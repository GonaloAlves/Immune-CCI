import os
import gc
import pandas as pd
import scanpy as sc
import gseapy as gp
import glob
import multiprocessing as mp

# Directory paths for t-statistic files and GSEA output
tstat_files = '/home/makowlg/Documents/Immune-CCI/src/tstat_files' 
gsea_dir = '/home/makowlg/Documents/Immune-CCI/src/gsea_dir'  # Directory for GSEA results

# Ensure the output directory exists
if not os.path.exists(gsea_dir):
    os.makedirs(gsea_dir)

# List all t-statistic Excel files in the directory
tstat_files = glob.glob(os.path.join(tstat_files, "*.xlsx"))

# Define gene set libraries to use (default for GSEA)
gene_sets = {
    "H:hallmark": ["MSigDB_Hallmark_2020"],
    "C2:curarted": ["KEGG_2016"],
    "C5:ontology": ["GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"],
} 

# Function to run GSEA using t-statistic rankings for a specific cluster


def prerank_gsea(kwargs: dict) -> None:
    """
    Perform GSEA analysis using the t-statistics for a specific cluster.
    The function loads gene names and their scores (t-statistics), performs GSEA, and outputs results.

    Parameters:
    kwargs (dict): Contains the following keys:
        - 'dataset': The dataset name.
        - 'cluster': The cluster name.
        - 'scores_df': The DataFrame containing the gene names and t-statistics.
        - 'gsea_dir': Directory for GSEA results.
        - 'threads': Number of threads for parallel processing.
        - 'gene_sets': Dictionary of gene sets for GSEA.

    Returns:
    None
    """

    seed = 2015854237  # Random seed for reproducibility
    dataset = kwargs['dataset']
    cluster = kwargs['cluster']
    
    print(f"Calculating GSEA for dataset {dataset} and cluster {cluster}...")

    # Sort by t-statistic values (scores)
    kwargs['scores_df'].sort_values(by=f'{cluster}.scores', ascending=False, inplace=True)
    kwargs['scores_df'].set_index(f'{cluster}.names', verify_integrity=True, inplace=True)

    # Build prerank DataFrame directly from 'scores_df'
    prerank_df = pd.DataFrame()
    prerank_df['gene'] = kwargs['scores_df'].index  # Directly take gene names from 'scores_df' index
    prerank_df['score'] = kwargs['scores_df'][f'{cluster}.scores']  # Take the t-statistic scores
    prerank_df.set_index('gene', verify_integrity=True, inplace=True)

    # Run GSEA for each gene set
    for gene_set_name, gene_set_paths in kwargs['gene_sets'].items():
        for gene_set_path in gene_set_paths:
            print(f"Running GSEA for {dataset}, cluster {cluster}, gene set {gene_set_name}...")

            # Define the output directory for the GSEA results
            gsea_output_dir = f"{kwargs['gsea_dir']}/{dataset}_{cluster}_{gene_set_name}"

            # Run GSEA
            gp.prerank(
                rnk=prerank_df,
                gene_sets=gene_set_path,
                threads=kwargs['threads'],
                min_size=10,   # Minimum size of gene sets
                max_size=2000,  # Maximum size of gene sets
                permutation_num=1000,  # Number of permutations
                graph_num=50,   # Number of graphs to show
                outdir=gsea_output_dir,  # Output directory
                seed=seed,
                verbose=True
            )

            # Clean up memory after each run
            gc.collect()

    print(f"GSEA for {dataset} - cluster {cluster} complete!")

def gsea_enrichment(dataset: str,
                    selected_genes: pd.DataFrame,
                    neighbor: int,
                    resolution: float,
                    enrichment_dir: str,
                    gsea_dir: str,
                    n_proc: int) -> None:
    dest = f"{enrichment_dir}/enrichment_mtx_{dataset}_final_leiden_n{neighbor}_r{resolution}.txt"
    print("Load gene enrichment data...")
    print(dest)
    df_enrich = pd.read_csv(dest, sep='\t', index_col=0, header=0)
    # Find clusters
    clusters = df_enrich.columns.to_list()
    clusters.sort(key=lambda x: int(x))
    
    # Run GSEA for each cluster
    print("Run pre-ranked GSEA...")
    n_tasks = len(clusters)
    print(f"Tasks to execute: {n_tasks}")
    if n_proc is None:
        n_proc = os.cpu_count()
    for c in clusters:
        args = {'scores_df': pd.DataFrame({f'{c}.names': df_enrich.loc[:, c].index.to_list(),
                                           f'{c}.scores': df_enrich.loc[:, c].values.tolist()}),
                'converted_genes_df': selected_genes,
                'dataset': dataset,
                'neighbor': neighbor,
                'resolution': resolution,
                'cluster': c,
                'gsea_dir': f'{gsea_dir}/enrichment_',
                'threads': n_proc,
                }
        prerank_gsea(args)
        
    # Free memory
    gc.collect()
    

def gsea_tstat(dataset: str,
               selected_genes: pd.DataFrame,
               neighbor: int,
               resolution: float,
               clusters: list,
               gsea_dir: str,
               n_proc: int) -> None:
    
    print("Load tstat tables...")
    tstat_dfs = {}
    for c in clusters:
        dest = f"{gsea_dir}/{dataset}_leiden_n{neighbor}_r{resolution}_c{c}_tstat.txt"
        tstat_dfs[c] = pd.read_csv(dest, sep='\t', index_col=0, header=0)
        print(dest)
        print(tstat_dfs[c])
    
    # Run GSEA for each cluster
    print("Run pre-ranked GSEA...")
    n_tasks = len(clusters)
    print(f"Tasks to execute: {n_tasks}")
    if n_proc is None:
        n_proc = os.cpu_count()
    for c in clusters:
        args = {'scores_df': pd.DataFrame({f'{c}.names': tstat_dfs[c].index.to_list(),
                                           f'{c}.scores': tstat_dfs[c]['t'].values.tolist()}),
                'converted_genes_df': selected_genes,
                'dataset': dataset,
                'neighbor': neighbor,
                'resolution': resolution,
                'cluster': c,
                'gsea_dir': f'{gsea_dir}/tstat_',
                'threads': n_proc,
                }
        prerank_gsea(args)
        
    # Free memory
    gc.collect()




# def run_gsea_for_cluster(tstat_file, gene_sets, output_dir):
#     """
#     Run GSEA analysis for a given t-statistic file.

#     Parameters:
#     - tstat_file (str): Path to the t-stat Excel file.
#     - gene_sets (str): Gene set database to use in the GSEA analysis.
#     - output_dir (str): Directory to save the GSEA results.
#     """
#     # Load the t-stat file into a DataFrame
#     df = pd.read_excel(tstat_file)

#     # Prepare the data for GSEA analysis
#     # Extract Gene Names and t-statistics for ranking
#     gsea_data = df[['Gene.Name', 't']]  # Replace 't' with the correct column containing t-statistics
    
#     # Sort by t-statistic for ranking
#     gsea_data = gsea_data.sort_values(by='t', ascending=False)

#     # Define output directory for the cluster
#     cluster_name = os.path.basename(tstat_file).split("_tstat")[0]
#     cluster_output_dir = os.path.join(output_dir, cluster_name)
#     if not os.path.exists(cluster_output_dir):
#         os.makedirs(cluster_output_dir)

#     # Run GSEA analysis using the t-statistics
#     print(f"Running GSEA for cluster {cluster_name}...")
#     results = gp.prerank(
#         rnk=gsea_data.values,  # Use Gene Name and t-statistics as ranking
#         gene_sets=gene_sets,
#         outdir=cluster_output_dir,  # Directory to save GSEA results
#         min_size=15,        # Minimum gene set size to consider
#         max_size=500,       # Maximum gene set size to consider
#         permutation_num=100,  # Number of permutations
#         seed=42  # For reproducibility
#     )
    
#     print(f"GSEA results saved for {tstat_file} in {cluster_output_dir}")
#     # Clean up memory
#     gc.collect()


# # Main function to process all t-stat files
# def start_gsea_tstat():
#     """
#     Start GSEA analysis for all clusters using available t-statistic files.
#     """
#     # Iterate through each t-statistic file and run GSEA analysis
#     for tstat_file in tstat_files:
#         # Run GSEA for each t-stat file
#         run_gsea_for_cluster(tstat_file, gene_sets, gsea_dir)


def start(n_proc=None) -> None:
    import mygene
    import time
    import datetime

    from src.globals import checkpoint_dir, annotation_dir, marker_genes_dir, gsea_dir
    from src.globals import data, datasets_divided, lineage_resolution_final, n_neighbors_final, n_hvg_spi
    
    neigh = n_neighbors_final[0]
    # Load socres
    t1 = time.time()
    for d in datasets_divided:
        dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked.h5ad"
        if os.path.exists(dest):
            print("Load gene rank data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)
        else:
            continue
        
        # Convert mouse genes to humam
        t2 = time.time()
        mg = mygene.MyGeneInfo()
        converted = mg.querymany(data[d].var_names, scopes='symbol,alias', species='human', fields='ensembl.gene,symbol', as_dataframe=True)
        converted.dropna(axis=0, subset='symbol', inplace=True)
        converted = converted['symbol']
        
        
        # Moderate t-statistic
        print("Calculate GSEA using R limma's moderate t-statistic...")
        clusters = data[d].obs[f'leiden_n_r{lineage_resolution_final[d][0]}'].cat.categories.to_list()
        gsea_tstat(dataset=d,
                  selected_genes=converted,
                  neighbor=neigh,
                  resolution=lineage_resolution_final[d][0],
                  clusters=clusters,
                  gsea_dir=gsea_dir,
                  n_proc=n_proc)
       
        
        print(f"Time for {d}: {datetime.timedelta(seconds=(time.time()-t2))}")
        
    print(f"Time total: {datetime.timedelta(seconds=(time.time()-t1))}")


# Main guard for multiprocessing support
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=os.cpu_count())

        print("\n********\n* DONE *\n********")

    except RuntimeError as e:
        print(f"An error occurred: {e}")