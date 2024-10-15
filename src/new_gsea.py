import os
import gc  # Garbage collector to free up memory
import pandas as pd
import scanpy as sc 
import gseapy as gp  # GSEApy is a Python wrapper for GSEA (Gene Set Enrichment Analysis)
import multiprocessing as mp  # For parallel processing

# Directory paths for t-statistic files and GSEA output
gsea_dir = '/home/makowlg/Documents/Immune-CCI/src/gsea_dir/Immune'  # Directory for GSEA results
tstat_dir = '/home/makowlg/Documents/Immune-CCI/src/tstat_files/Immune'  # Directory for t-statistic files

# Ensure the output directory exists
if not os.path.exists(gsea_dir):
    os.makedirs(gsea_dir)  # Create the directory if it doesn't exist

# Define gene set libraries to use (default for GSEA)
gene_sets = {
    "H:hallmark": ["MSigDB_Hallmark_2020"],  # Hallmark pathways
    "C2:curated": ["KEGG_2016"],  # KEGG curated pathways
    "C5:ontology": ["GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"],  # Gene Ontology
}

# Function to run GSEA using t-statistic rankings for a specific cluster
def prerank_gsea(kwargs: dict) -> None:
    """
    Run pre-ranked GSEA for a specific cluster.

    Parameters:
    kwargs (dict): A dictionary containing the parameters for GSEA.
    """
    
    seed = 2015854237  # Seet seed for reproducibility
    dataset = kwargs['dataset']  # The dataset name
    cluster = kwargs['cluster']  # The specific cluster to process
    resolution = kwargs['resolution']  # The resolution used for clustering
    
    print(f"Calculating GSEA for dataset {dataset} and cluster {cluster}...")

    # Sort by t-statistic values (scores)
    kwargs['scores_df'].sort_values(by=f'{cluster}.scores', ascending=False, inplace=True)
    kwargs['scores_df'].set_index(f'{cluster}.names', verify_integrity=True, inplace=True)  # Set gene names as index

    # Build a DataFrame for prerank input (gene and scores)
    prerank_df = pd.DataFrame()
    prerank_df['gene'] = kwargs['converted_genes_df'].loc[kwargs['scores_df'].index[kwargs['scores_df'].index.isin(kwargs['converted_genes_df'].index)]]
    prerank_df.drop_duplicates(subset='gene', inplace=True)  # Remove duplicate genes
    prerank_df['score'] = kwargs['scores_df'].loc[prerank_df.index, f'{cluster}.scores']
    prerank_df.set_index('gene', verify_integrity=True, inplace=True)  # Set gene as index

    # Run GSEA for each gene set defined earlier
    for set_collection, set_collection_name in gene_sets.items():
        for gene_set_path in set_collection_name:
            print(f"Analyzing {dataset} leiden_{resolution}_c{cluster} {set_collection} - {gene_set_path}...")
            # Set destination path for GSEA results
            dest = kwargs['gsea_dir'] + set_collection.replace(':', '_') + f"/{dataset}_leiden_{resolution}_c{cluster}_{gene_set_path}"

            # Perform pre-ranked GSEA using the prerank DataFrame
            pre_res = gp.prerank(
                rnk=prerank_df,
                gene_sets=gene_set_path,
                threads=kwargs['threads'],
                min_size=10,  # Minimum gene set size for enrichment analysis
                max_size=2000,  # Maximum gene set size
                permutation_num=1000,  # Number of permutations for significance testing
                graph_num=50,  # Number of graphs to generate
                outdir=dest,  # Output directory for results    
                seed=seed,  # Set seed for reproducibility
                verbose=True  # Print detailed output
            )

            # Free memory after each GSEA run
            gc.collect()

    print(f"GSEA for {dataset} - cluster {cluster} complete!")

# Function to load t-statistic data and run GSEA for each cluster
def gsea_tstat(dataset: str,
               selected_genes: pd.DataFrame,
               resolution: str,
               clusters: list,
               gsea_dir: str,
               tstat_dir: str,
               n_proc: int) -> None:
    """
    Load t-statistic data and run pre-ranked GSEA for each cluster.

    Parameters:
    dataset (str): Name of the dataset.
    selected_genes (pd.DataFrame): DataFrame containing selected genes.
    resolution (float): Clustering resolution.
    clusters (list): List of clusters to process.
    gsea_dir (str): Path to the directory for GSEA output.
    tstat_dir (str): Path to the directory containing t-statistic files.
    n_proc (int): Number of parallel processes for multiprocessing.
    """
    
    print("Load tstat tables...")
    tstat_dfs = {}  # Dictionary to store t-statistic DataFrames
    
    # Iterate over clusters and load the corresponding t-statistic data
    for c in clusters:
        dest = f"{tstat_dir}/{dataset}_leiden_{resolution}_c{c}_tstat.xlsx"
        tstat_dfs[c] = pd.read_csv(dest, sep='\t', index_col=0, header=0)  # Load t-statistic data
        print(tstat_dfs[c])  # Print the loaded data for inspection
    
    # Run GSEA for each cluster
    print("Run pre-ranked GSEA...")
    n_tasks = len(clusters)
    print(f"Tasks to execute: {n_tasks}")
    if n_proc is None:  # If n_proc is not set, use all available CPU cores
        n_proc = os.cpu_count()

    # Loop over each cluster and run GSEA
    for c in clusters:
        args = {
            'scores_df': pd.DataFrame({f'{c}.names': tstat_dfs[c].index.to_list(),
                                       f'{c}.scores': tstat_dfs[c]['t'].values.tolist()}),  # Prepare DataFrame with names and scores
            'converted_genes_df': selected_genes,  # DataFrame with converted genes
            'dataset': dataset,  # Dataset name
            'resolution': resolution,  # Clustering resolution
            'cluster': c,  # Cluster to process
            'gsea_dir': f'{gsea_dir}/tstat_',  # Directory for GSEA output
            'threads': n_proc,  # Number of threads (processes) for GSEA
        }
        prerank_gsea(args)  # Call the GSEA function
    
    # Free memory after processing
    gc.collect()

# Entry function for the entire workflow
def start(n_proc=None) -> None:
    """
    Start the GSEA analysis pipeline, including loading data, converting genes, and running GSEA.

    Parameters:
    n_proc (int, optional): Number of processes to use for multiprocessing.
    """
    import mygene  # Module for gene conversion (mouse to human)
    import time  # Time tracking for performance measurement
    import datetime  # Time formatting for output
    
    # Load scores from an h5ad file
    t1 = time.time()  # Start time tracking
    
    dest = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad"
    if os.path.exists(dest):  # Check if the file exists
        print("Load gene rank data...")
        data = sc.read_h5ad(dest)  # Load the AnnData object containing the data
        
        
        # Convert mouse genes to human genes using mygene
        t2 = time.time()
        mg = mygene.MyGeneInfo()  # Initialize mygene API
        converted = mg.querymany(data.var_names, # list of gene names 
                                 scopes='symbol,alias', # look for gene symbols and aliases ensuring that all recognized names are considered.
                                 species='human', # Return human gene information
                                 fields='ensembl.gene,symbol', # Ensembl gene ID for the human gene | Official gene symbol
                                 as_dataframe=True # Df return
        )  
        converted.dropna(axis=0, 
                         subset='symbol', 
                         inplace=True)  # Remove rows with missing gene symbols
        
        converted = converted['symbol']  # Select only the 'symbol' column
        
        # Perform GSEA using the moderate t-statistic
        print("Calculate GSEA using R limma's moderate t-statistic...")
        
        dataset = 'Immune'  # Set the dataset name
        clusters = data.obs['leiden_fusion'].cat.categories.to_list()  # Get the list of clusters
        
        gsea_tstat(dataset=dataset,
                   selected_genes=converted,
                   resolution='fusion',  # Use 'fusion' as resolution
                   clusters=clusters,  # List of clusters
                   gsea_dir=gsea_dir,
                   tstat_dir=tstat_dir,
                   n_proc=n_proc)  # Number of processes to use
        
        print(f"Time for {dataset}: {datetime.timedelta(seconds=(time.time()-t2))}")  # Print time for GSEA calculation
    
    print(f"Time total: {datetime.timedelta(seconds=(time.time()-t1))}")  # Print total execution time

# Main guard to enable multiprocessing (especially for Windows)
if __name__ == '__main__':
    try:
        mp.set_start_method('spawn', force=True)  # Ensure compatibility with Windows for multiprocessing
        start(n_proc=os.cpu_count())  # Call the start function with the available CPU core count

        print("\n********\n* DONE *\n********")  # Indicate completion

    except RuntimeError as e:
        print(f"An error occurred: {e}")  # Catch and print any runtime errors
