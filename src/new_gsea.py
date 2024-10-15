import os
import gc
import pandas as pd
import scanpy as sc
import gseapy as gp
import multiprocessing as mp

# Directory paths for t-statistic files and GSEA output
gsea_dir = '/home/makowlg/Documents/Immune-CCI/src/gsea_dir/Immune'  # Directory for GSEA results
tstat_dir = '/home/makowlg/Documents/Immune-CCI/src/tstat_files/Immune'

# Ensure the output directory exists
if not os.path.exists(gsea_dir):
    os.makedirs(gsea_dir)

# Define gene set libraries to use (default for GSEA)
gene_sets = {
    "H:hallmark": ["MSigDB_Hallmark_2020"],
    "C2:curarted": ["KEGG_2016"],
    "C5:ontology": ["GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"],
} 

# Function to run GSEA using t-statistic rankings for a specific cluster
def prerank_gsea(kwargs: dict) -> None:

    seed = 2015854237  # Random seed for reproducibility
    dataset = kwargs['dataset']
    cluster = kwargs['cluster']
    resolution = kwargs['resolution']
    
    
    print(f"Calculating GSEA for dataset {dataset} and cluster {cluster}...")

    # Sort by t-statistic values (scores)
    kwargs['scores_df'].sort_values(by=f'{cluster}.scores', ascending=False, inplace=True)
    kwargs['scores_df'].set_index(f'{cluster}.names', verify_integrity=True, inplace=True)

    # Build prerank df
    prerank_df = pd.DataFrame()
    prerank_df['gene'] = kwargs['converted_genes_df'].loc[kwargs['scores_df'].index[kwargs['scores_df'].index.isin(kwargs['converted_genes_df'].index)]]
    prerank_df.drop_duplicates(subset='gene', inplace=True)
    prerank_df['score'] = kwargs['scores_df'].loc[prerank_df.index, f'{cluster}.scores']
    prerank_df.set_index('gene', verify_integrity=True, inplace=True)
        

    # Run GSEA for each gene set
    for set_collection, set_collection_name in gene_sets.items():
        for gene_set_path in set_collection_name:
            print(f"Analyzing {dataset} leiden_{resolution}_c{cluster} {set_collection} - {gene_set_path}...")
            dest = kwargs['gsea_dir'] + set_collection.replace(':', '_') + f"/{dataset}_leiden_{resolution}_c{cluster}_{gene_set_path}"
            pre_res = gp.prerank(rnk=prerank_df,
                                    gene_sets=gene_set_path,
                                    threads=kwargs['threads'],
                                    min_size=10,   # 5
                                    max_size=2000,  # [500-1000]
                                    permutation_num=1000,  # reduce number to speed up testing
                                    graph_num=50,   # 20
                                    outdir=dest,
                                    seed=seed,
                                    verbose=True)
            #print(pre_res)
            
            # Free memory
            gc.collect()

    print(f"GSEA for {dataset} - cluster {cluster} complete!")


def gsea_tstat(dataset: str,
               selected_genes: pd.DataFrame,
               resolution: float,
               clusters: list,
               gsea_dir: str,
               tstat_dir: str,
               n_proc: int) -> None:
    
    print("Load tstat tables...")
    tstat_dfs = {}
    
    print(dataset)
    for c in clusters:
        dest = f"{tstat_dir}/{dataset}_leiden_{resolution}_c{c}_tstat.xlsx"
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
                'resolution': resolution,
                'cluster': c,
                'gsea_dir': f'{gsea_dir}/tstat_',
                'threads': n_proc,
                }
        prerank_gsea(args)
        
    # Free memory
    gc.collect()



def start(n_proc=None) -> None:
    import mygene
    import time
    import datetime

    # Load socres
    t1 = time.time()
    
    dest = "/home/makowlg/Documents/Immune-CCI/h5ad_files/adata_final_Immune_raw_norm_ranked_copy.h5ad"
    if os.path.exists(dest):
        print("Load gene rank data...")
        data = sc.read_h5ad(dest)
        dataset = 'Immune'
        
        # Convert mouse genes to humam
        t2 = time.time()
        mg = mygene.MyGeneInfo()
        converted = mg.querymany(data.var_names, scopes='symbol,alias', species='human', fields='ensembl.gene,symbol', as_dataframe=True)
        converted.dropna(axis=0, subset='symbol', inplace=True)
        converted = converted['symbol']
        
        
        # Moderate t-statistic
        print("Calculate GSEA using R limma's moderate t-statistic...")
        clusters = data.obs['leiden_fusion'].cat.categories.to_list()
        
        gsea_tstat(dataset=dataset,
                  selected_genes=converted,
                  resolution='fusion',
                  clusters=clusters,
                  gsea_dir=gsea_dir,
                  tstat_dir = tstat_dir,
                  n_proc=n_proc)
       
        
        print(f"Time for {dataset}: {datetime.timedelta(seconds=(time.time()-t2))}")
        
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