# Gene Set enrichment analysis (GSEA) - LeonorSaude10x
#
# Build a summary of the results from S7
#
# Daniel Ribeiro, 2023
import os

import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('cairo')


def collect_summary(df: pd.DataFrame,
                    pval: float = 0.05) -> pd.DataFrame:
    result = df[df['FDR q-val'] < pval]
    result = result.loc[:, ['Term', 'NES', 'FDR q-val']]
    result.sort_values(by='NES', ascending=False, inplace=True)
    
    return result
    

def collection_summary(dataset: str,
                       geneset: str,
                       clusters: list,
                       collection_path: str) -> pd.DataFrame:
    summary = pd.DataFrame()
    print("Load GSEA data...")
    for c in clusters:
        dest = f"{collection_path}/{dataset}_weighted_SC_{c}_{geneset}/gseapy.gene_set.prerank.report.csv"      # TODO: change name
        if os.path.exists(dest):
            df = pd.read_csv(dest, sep=',', header=0, index_col=None)
            print(dest)
        else:
            continue
        
        df = collect_summary(df)
        df.reset_index(inplace=True)
        df['index'] = df.index.to_list()
        if not df.empty:
            # Get cluster number
            cluster_n = c.split('c')[-1]
            cols = df.columns
            cols = [f'{cluster_n}.{col}' for col in cols]
            df.columns = cols
            summary = pd.concat([summary, df], axis=1)
    summary.fillna('', inplace=True)
    
    return summary


def start() -> None:
    import src.globals   # noqa:F401
    from src.globals import checkpoint_dir, gsea_dir
    from src.globals import data, datasets_divided, senescence_resolution, n_neighbors_final, gene_sets
    
    neigh = n_neighbors_final[0]
    # Load scores
    for d in datasets_divided:
        dest = f"{checkpoint_dir}/adata_final_{d}_weighted_SC_raw_norm_ranked_genes.h5ad"       # TODO: change name
        if os.path.exists(dest):
            print("Load gene rank data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)
        else:
            continue
        
        # Find cluster names
        if d.find('double') != -1:
            integration = 'NES_INTEGRATED_DOUBLE'
        elif d.find('triple') != -1:
            integration = 'NES_INTEGRATED_TRIPLE'
        else:
            integration = 'NES_INTEGRATED'
        cluster = f"leiden_n{neigh}_r{senescence_resolution[dataset][integration][0]}"      # TODO: get correct cluster resolution
        cluster_list = data[d].obs[cluster].cat.categories.tolist()
        cluster_list = [f'{cluster}_c{c}' for c in cluster_list]
        
        # Collect summary analysis
        scores = ["wilcoxon", "tstat"]           # TODO: simplify, since only tstat was used
        for score in scores:
            for collection, gsets in gene_sets.items():
                dest = f"{gsea_dir}/{score}_{collection.replace(':', '_')}"
                for g in gsets:
                    summary = collection_summary(dataset=d,
                                                 geneset=g,
                                                 clusters=cluster_list,
                                                 collection_path=dest)
                    if not summary.empty:
                        print("Save summary data...")
                        dest2 = f"{gsea_dir}/summary_{score}_{d}_weighted_SC_{cluster}_{g}.txt"
                        summary.to_csv(dest2, sep='\t', index=False)
                        print(dest2)
                        print()
            
    
start()

print("\n********\n* DONE *\n********")
