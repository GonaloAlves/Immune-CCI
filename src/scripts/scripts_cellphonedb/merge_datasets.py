# Merge the divided datsets into 1 - LeonorSaude10x
#
# Also, rank genes
#
# Daniel Ribeiro, 2023
import gc
import os

import scanpy as sc
import pandas as pd
    

cellphonedb = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb"
checkpoint_dir = "/home/makowlg/Documents/Immune-CCI/h5ad_files"

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
    marker_genes.to_csv(f'{cellphonedb}/{prefix}_{name}.csv', sep='\t')


def start() -> None:

    data: dict[str, sc.AnnData] = {}
    datasets_divided = [
    'Neu_CentralCanal',
    'Meningeal_Vascular',
    'Immune'
    ]
    # .obs to merge
    obs_merge = ['injury', 'day', 'collection_region', 'injury_day', 'injury_region', 'injury_condition',
                 'nuclei_uL', 'total_nuclei', 'target_10x', 'mouse_id', 'date_nuclei_extraction', 'sample_id',
                 'sample_I15AC', 'sample_I15AR', 'sample_I15BC', 'sample_I15CC', 'sample_I15CR', 'sample_I60AC',
                 'sample_I60AR', 'sample_I60BC', 'sample_I60BR', 'sample_I60CC', 'sample_I60CR', 'sample_S15AC',
                 'sample_S15BC', 'sample_S15BR', 'sample_S15CC', 'sample_S15CR', 'sample_U00AX', 'sample_U00BX',
                 'sample_U00CX', 'n_counts', 'filt_counts', 'n_genes', 'filt_genes', 'percent_mito', 'filt_mito',
                 'doublet_score', 'predicted_doublet', 'doublet', 'annot_lineage_cell']
    
    datset_names = {'Neu_CentralCanal': 'Neu',
                    'Meningeal_Vascular': 'MeV',
                    'Immune': 'Imm'
                    }
    obs_list = []
    
    # Load scores
    for d in datasets_divided:
        dest = f"{checkpoint_dir}/adata_final_{d}_raw_norm_ranked_copy_copy.h5ad"
        if os.path.exists(dest):
            print("Load gene rank data...")
            print(dest)
            data[d] = sc.read_h5ad(dest)
            data[d].uns['log1p']['base'] = None
            
            # Name of .obs column to access
            obs_list.append(f"leiden_fusion")
            
            # Rename categories
            data[d].obs['leiden_merge'] = data[d].obs[obs_list[-1]]
            # data[d].obs['leiden_merge'] = data[d].obs['leiden_merge'].cat.rename_categories(lambda x: f'{datset_names[d]}.{x}')
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
                    
    
    # Filter rows where 'leiden_merge' ends with 'NA'
    adata_concat = adata_concat[~adata_concat.obs['leiden_merge'].str.endswith('NA')]
    
    
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
    print("All cells--------------------------")
    adata_concat.uns['rank_genes_groups_leiden_merge'] = result
    write_marker_genes(adata_concat, rank_genes_groups='rank_genes_groups_leiden_merge', prefix="adata_final_merged_marker_genes")
    print("END All cells--------------------------")
    print(adata_concat.obs['leiden_merge'].value_counts())

    # for subset in ["uninjured.0", "sham.15", "injured.15", "injured.60"]:
    #     subset_data = adata_concat[adata_concat.obs['injury_day'] == subset]
    #     print(f"\nSubset: {subset}, Total Cells: {subset_data.shape[0]}")
    #     print(subset_data.obs['leiden_merge'].value_counts())

    # print("END All cells--------------------------")
    # # Uninjured cells
    # subset = "uninjured.0"
    # args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
    #        'n_genes': None,
    #        'groupby': 'leiden_merge',
    #        'method': 'wilcoxon',
    #        'use_raw': False,
    #        'key_added': f'rank_genes_groups_leiden_merge_{subset}',
    #        'pts': True}
    # result = rank_genes(args)
    # print("Uninjured --------------------------")
    # print(result)
    # print("Uninjured --------------------------")
    # adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    # write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")

    # # Sham cells
    # subset = "sham.15"
    # args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
    #         'n_genes': None,
    #         'groupby': 'leiden_merge',
    #         'method': 'wilcoxon',
    #         'use_raw': False,
    #         'key_added': f'rank_genes_groups_leiden_merge_{subset}',
    #         'pts': True}
    # result = rank_genes(args)
    # print("sham 15--------------------------")
    # print(result)
    # print("sham 15--------------------------")
    # adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    # write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")

    # # Injury_15 cells
    # subset = "injured.15"
    # args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
    #         'n_genes': None,
    #         'groupby': 'leiden_merge',
    #         'method': 'wilcoxon',
    #         'use_raw': False,
    #         'key_added': f'rank_genes_groups_leiden_merge_{subset}',
    #         'pts': True}
    # result = rank_genes(args)
    # print("Injured 15--------------------------")
    # print(result)
    # print("Injured 15--------------------------")
    # adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    # write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")

    # # Injury_60 cells
    # subset = "injured.60"
    # args = {'adata': adata_concat[adata_concat.obs['injury_day'] == subset].copy(),
    #         'n_genes': None,
    #         'groupby': 'leiden_merge',
    #         'method': 'wilcoxon',
    #         'use_raw': False,
    #         'key_added': f'rank_genes_groups_leiden_merge_{subset}',
    #         'pts': True}
    # result = rank_genes(args)
    # print("Injured 60--------------------------")
    # print(result)
    # print("Injured 60--------------------------")
    # adata_concat.uns[f'rank_genes_groups_leiden_merge_{subset}'] = result
    # write_marker_genes(adata_concat, rank_genes_groups=f'rank_genes_groups_leiden_merge_{subset}', prefix="adata_final_merged_marker_genes")
    
    
    # Save AnnData
    print("Save merged AnnData...")
    dest = f"{checkpoint_dir}/adata_final_merged_raw_norm_annot_nona.h5ad"
    print(dest)
    adata_concat.write_h5ad(dest, compression='gzip')
        

start()

print("\n********\n* DONE *\n********")
