# Merge the divided datsets into 1 - LeonorSaude10x
#
# Also, rank genes
#
# Daniel Ribeiro, 2023
import gc
import os

import scanpy as sc
import pandas as pd
    

def rank_genes(kwargs: dict):
    print(f"Ranking genes on {kwargs['key_added']}...")
    if len(kwargs['adata'].obs[kwargs['groupby']].cat.categories.tolist()) > 1:
        sc.tl.rank_genes_groups(**kwargs)
        gc.collect()
        return kwargs['adata'].uns[kwargs['key_added']].copy()
    else:
        return None


def write_marker_genes(adata: sc.AnnData,
                       rank_genes_groups: str,      # Key in .uns
                       prefix: str                  # File name prefix
                       ) -> None:
    from src.globals import marker_genes_dir
    
    print("Write marker genes...")
    marker_genes = pd.DataFrame()
    # Get size col numbers
    n_col = len(pd.DataFrame(adata.uns[rank_genes_groups]['names']).columns)
    categories = pd.DataFrame(adata.uns[rank_genes_groups]['names']).columns.to_list()
    n_row = len(pd.DataFrame(adata.uns[rank_genes_groups]['names']).index)
    df = {}
    for col in range(n_col):
        # Categories are ordered alphabetically. Order is low, top
        df_temp = pd.DataFrame()
        df_temp[f'{categories[col]}.names'] = pd.DataFrame(adata.uns[rank_genes_groups]['names']).iloc[:n_row, col]
        df_temp[f'{categories[col]}.scores'] = pd.DataFrame(adata.uns[rank_genes_groups]['scores']).iloc[:n_row, col]
        df_temp[f'{categories[col]}.logfoldchanges'] = pd.DataFrame(adata.uns[rank_genes_groups]['logfoldchanges']).iloc[:n_row, col]
        df_temp[f'{categories[col]}.pvals'] = pd.DataFrame(adata.uns[rank_genes_groups]['pvals']).iloc[:n_row, col]
        df_temp[f'{categories[col]}.pvals_adj'] = pd.DataFrame(adata.uns[rank_genes_groups]['pvals_adj']).iloc[:n_row, col]
        df_temp[f'{categories[col]}.pts'] = pd.DataFrame(adata.uns[rank_genes_groups]['pts']).loc[df_temp[f'{categories[col]}.names'], categories[col]].values
        df_temp.sort_values(by=f'{categories[col]}.logfoldchanges', ascending=False, ignore_index=True, inplace=True)
        df_temp = df_temp.copy()
        
        df[f'{categories[col]}.names'] = df_temp[f'{categories[col]}.names']
        df[f'{categories[col]}.scores'] = df_temp[f'{categories[col]}.scores']
        df[f'{categories[col]}.logfoldchanges'] = df_temp[f'{categories[col]}.logfoldchanges']
        df[f'{categories[col]}.pvals'] = df_temp[f'{categories[col]}.pvals']
        df[f'{categories[col]}.pvals_adj'] = df_temp[f'{categories[col]}.pvals_adj']
        df[f'{categories[col]}.pts'] = df_temp[f'{categories[col]}.pts']
    
    # Save table
    name = rank_genes_groups[18:]   # The first 18chars are always the same
    marker_genes = pd.DataFrame(df)
    marker_genes = marker_genes.copy()
    marker_genes.to_csv(f'{marker_genes_dir}/{prefix}_{name}.txt', sep='\t')


def start() -> None:
    import src.globals   # noqa:F401
    from src.globals import checkpoint_dir
    from src.globals import data, datasets_divided, n_neighbors_final, lineage_resolution_final
    
    # .obs to merge
    obs_merge = ['injury', 'day', 'collection_region', 'injury_day', 'injury_region', 'injury_condition',
                 'nuclei_uL', 'total_nuclei', 'target_10x', 'mouse_id', 'date_nuclei_extraction', 'sample_id',
                 'sample_I15AC', 'sample_I15AR', 'sample_I15BC', 'sample_I15CC', 'sample_I15CR', 'sample_I60AC',
                 'sample_I60AR', 'sample_I60BC', 'sample_I60BR', 'sample_I60CC', 'sample_I60CR', 'sample_S15AC',
                 'sample_S15BC', 'sample_S15BR', 'sample_S15CC', 'sample_S15CR', 'sample_U00AX', 'sample_U00BX',
                 'sample_U00CX', 'n_counts', 'filt_counts', 'n_genes', 'filt_genes', 'percent_mito', 'filt_mito',
                 'doublet_score', 'predicted_doublet', 'doublet', 'annot_lineage_cell']
    datset_names = {'Neuron': 'Neu',
                    'Oligodendrocyte': 'Oli',
                    'Astrocyte': 'Ast',
                    'Meningeal_Vascular': 'MeV',
                    'Immune': 'Imm',
                    }
    obs_list = []
    
    # Load scores
    for d in datasets_divided:
        dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked.h5ad"
        if os.path.exists(dest):
            print("Load gene rank data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)
            data[d].uns['log1p']['base'] = None
            
            # Name of .obs column to access
            obs_list.append(f"leiden_n{n_neighbors_final[0]}_r{lineage_resolution_final[d][0]}")
            
            # Rename categories
            data[d].obs['leiden_merge'] = data[d].obs[obs_list[-1]]
            data[d].obs['leiden_merge'] = data[d].obs['leiden_merge'].cat.rename_categories(lambda x: f'{datset_names[d]}.{x}')
            # scDiffComm cannot does not suppor '_' in category names --> substitute them
            data[d].obs['injury_day'] = data[d].obs['injury_day'].cat.rename_categories(lambda x: x.replace('_', '.'))
            data[d].obs['injury_region'] = data[d].obs['injury_region'].cat.rename_categories(lambda x: x.replace('_', '.'))
            data[d].obs['injury_condition'] = data[d].obs['injury_condition'].cat.rename_categories(lambda x: x.replace('_', '.'))
        else:
            continue
        
    adata_concat = sc.concat(adatas=data.values(), merge='same')
    # Only keep meaningful .obs
    obs_merge.append("leiden_merge")
    adata_keys = list(adata_concat.obs_keys())
    for k in adata_keys:
        if k not in obs_merge:
            del adata_concat.obs[k]
                    
    # Rank genes
    print("Rank genes...")
    # All cells
    args = {'adata': adata_concat,
            'n_genes': None,
            'groupby': 'leiden_merge',
            'method': 'wilcoxon',
            'use_raw': False,
            'key_added': 'rank_genes_groups_leiden_merge',
            'pts': True}
    result = rank_genes(args)
    adata_concat.uns['rank_genes_groups_leiden_merge'] = result
    write_marker_genes(adata_concat, rank_genes_groups='rank_genes_groups_leiden_merge', prefix="adata_final_merged_marker_genes")
    ## Uninjured cells
    #subset = "uninjured.0"
    #args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
    #        'n_genes': None,
    #        'groupby': 'leiden_merge',
    #        'method': 'wilcoxon',
    #        'use_raw': False,
    #        'key_added': f'rank_genes_groups_leiden_merge_{subset}',
    #        'pts': True}
    #result = rank_genes(args)
    #adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    #write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")
    # Sham cells
    subset = "sham.15"
    args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
            'n_genes': None,
            'groupby': 'leiden_merge',
            'method': 'wilcoxon',
            'use_raw': False,
            'key_added': f'rank_genes_groups_leiden_merge_{subset}',
            'pts': True}
    result = rank_genes(args)
    adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")
    # Injury_15 cells
    subset = "injured.15"
    args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
            'n_genes': None,
            'groupby': 'leiden_merge',
            'method': 'wilcoxon',
            'use_raw': False,
            'key_added': f'rank_genes_groups_leiden_merge_{subset}',
            'pts': True}
    result = rank_genes(args)
    adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")
    # Injury_60 cells
    subset = "injured.60"
    args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
            'n_genes': None,
            'groupby': 'leiden_merge',
            'method': 'wilcoxon',
            'use_raw': False,
            'key_added': f'rank_genes_groups_leiden_merge_{subset}',
            'pts': True}
    result = rank_genes(args)
    adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")
    
    ## Rostral-caudal cells
    #for d in datasets_divided:
    #    # collection_region
    #    args = {'adata': data[d],
    #            'n_genes': None,
    #            'groupby': 'collection_region',
    #            'method': 'wilcoxon',
    #            'use_raw': False,
    #            'key_added': 'rank_genes_groups_collection_region',
    #            'pts': True}
    #    result = rank_genes(args)
    #    data[d].uns['rank_genes_groups_collection_region'] = result
    #    write_marker_genes(data[d], rank_genes_groups=f'rank_genes_groups_collection_region_{subset}', prefix=f"adata_final_{d}_marker_genes")
    #    # injury_region
    #    args = {'adata': data[d],
    #            'n_genes': None,
    #            'groupby': 'injury_region',
    #            'method': 'wilcoxon',
    #            'use_raw': False,
    #            'key_added': 'rank_genes_groups_injury_region',
    #            'pts': True}
    #    result = rank_genes(args)
    #    data[d].uns['rank_genes_groups_injury_region'] = result
    #    write_marker_genes(data[d], rank_genes_groups=f'rank_genes_groups_injury_region_{subset}', prefix=f"adata_final_{d}_marker_genes")
    #    # injury_condition
    #    args = {'adata': data[d],
    #            'n_genes': None,
    #            'groupby': 'injury_condition',
    #            'method': 'wilcoxon',
    #            'use_raw': False,
    #            'key_added': 'rank_genes_groups_injury_condition',
    #            'pts': True}
    #    result = rank_genes(args)
    #    data[d].uns['rank_genes_groups_injury_condition'] = result
    #    write_marker_genes(data[d], rank_genes_groups=f'rank_genes_groups_injury_condition_{subset}', prefix=f"adata_final_{d}_marker_genes")
    #
    #    # Save AnnData
    #    print("Save AnnData...")
    #    dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked.h5ad"
    #    print(dest)
    #    adata_concat.write_h5ad(dest, compression='gzip')
    
    # Save AnnData
    print("Save merged AnnData...")
    dest = f"{checkpoint_dir}/adata_final_merged_raw_norm_annot.h5ad"
    print(dest)
    adata_concat.write_h5ad(dest, compression='gzip')
        

start()

print("\n********\n* DONE *\n********")
