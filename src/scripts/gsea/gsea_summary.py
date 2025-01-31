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

gsea_dir = '/home/makowlg/Documents/Immune-CCI/src/gsea_dir/Meningeal'  
h5ad_dir = '/home/makowlg/Documents/Immune-CCI/h5ad_files/'

dataset = 'Meningeal_Vascular'
resolution = 'fusion'

gene_sets = {
    "H:hallmark": ["MSigDB_Hallmark_2020"],  # Hallmark pathways
    "C2:curated": ["KEGG_2016"],  # KEGG curated pathways
    "C5:ontology": ["GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"],  # Gene Ontology
}


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
        dest = f"{collection_path}/{dataset}_{c}_{geneset}/gseapy.gene_set.prerank.report.csv"     
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
    # Load scores
    
    dest = f"{h5ad_dir}/adata_final_{dataset}_raw_norm_ranked_copy_copy.h5ad"    
    if os.path.exists(dest):
        print("Load gene rank data...")
        print(dest)
        data = sc.read_h5ad(dest)
    
    
    # Find cluster names
    if dataset.find('double') != -1:
        integration = 'NES_INTEGRATED_DOUBLE'
    elif dataset.find('triple') != -1:
        integration = 'NES_INTEGRATED_TRIPLE'
    else:
        integration = 'NES_INTEGRATED'
    cluster = f"leiden_{resolution}"      

    cluster_list = data.obs[cluster].cat.categories.tolist()
    cluster_list = [f'{cluster}_c{c}' for c in cluster_list]


    # Collect summary analysis
    scores = ["tstat"]           
    for score in scores:
        for collection, gsets in gene_sets.items():
            dest = f"{gsea_dir}/{score}_{collection.replace(':', '_')}"
            for g in gsets:
                summary = collection_summary(dataset=dataset,
                                            geneset=g,
                                            clusters=cluster_list,
                                            collection_path=dest)
                if not summary.empty:
                    print("Save summary data...")
                    dest2 = f"{gsea_dir}/summary_{score}_{dataset}_raw_norm_ranked_{cluster}_{g}.csv"
                    summary.to_csv(dest2, sep='\t', index=False)
                    print(dest2)
                    print()
        
    
start()

print("\n********\n* DONE *\n********")
