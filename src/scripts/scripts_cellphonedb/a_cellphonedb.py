# Perform CellPhoneDB analysis - LeonorSaude10x
#
# Follwing these tutorials
# 
# Daniel Ribeiro, Gonçalo Alves 2025
import pandas as pd
import scanpy as sc

simple_analysis = True
statistical_analysis = True
deg_analysis = True

cellphonedb = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb"
cpdb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database"
marker_genes_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/markers"
checkpoint_dir = "/home/makowlg/Documents/Immune-CCI/h5ad_files"

def build_meta_file(adata: sc.AnnData,
                    cell_type: str,        # .obs key containing cell_type data
                    group_by: str = None,  # .obs key to subset data. If None, use all cells
                    subset: str = None,    # subset key
                    ) -> str:
    
    if group_by is None:    # All cells
        print("Build metafile for all cells")
        meta_df = pd.DataFrame({"barcode_sample": adata.obs_names.to_list(),
                                "cell_type": adata.obs[cell_type].values.tolist()})
        meta_file_path = f"{cellphonedb}/adata_final_merged_metadata_nona.tsv"
        print(meta_file_path)
        meta_df.to_csv(meta_file_path, sep='\t', index=False)
        return meta_file_path
    elif (group_by is not None) and (subset is not None):
        print(f"Build metafile for {group_by} {subset} cells")
        subset_cells = adata[adata.obs[group_by] == subset].obs_names.to_list()
        meta_df = pd.DataFrame({"barcode_sample": adata[subset_cells, :].obs_names.to_list(),
                                "cell_type": adata[subset_cells, :].obs[cell_type].values.tolist()})
        meta_file_path = f"{cellphonedb}/adata_final_merged_{subset}_metadata_nona.tsv"
        print(meta_file_path)
        meta_df.to_csv(meta_file_path, sep='\t', index=False)
        return meta_file_path
    else:
        raise ValueError("'group_by' and 'subset' must be valid .obs keys and category, respectively.")


def build_deg_file(file_path: str,
                   out_file: str) -> str:
    import mygene
    
    # Load table
    print(f"\nBuild {out_file} deg file")
    df = pd.read_csv(file_path, header=0, index_col=0, sep='\t')
    # Separate cell types
    cell_types = []
    for col in df.columns:
        split = col.rsplit('.', maxsplit=1)[0]
        if split not in cell_types:
            cell_types.append(split)
    # Separate dfs of cell types
    df_separated = {}
    for cell in cell_types:
        mask = df.columns.str.startswith(cell)
        df_cell_type = df.loc[:, mask].copy()
        df_cell_type.sort_values(by=f'{cell}.names', inplace=True)    # Sort by gene
        df_cell_type.reset_index(drop=True, inplace=True)             # Reset index
        df_separated[cell] = df_cell_type
    
    # Convert ms genes to hs genes
    # The genes for all cell cell_types are the same
    # Use first cell type to get ms genes
    print("Convert mouse genes to human...")
    mg = mygene.MyGeneInfo()
    cell_type = cell_types[0]
    genes = df.loc[:, f'{cell_type}.names'].copy()
    genes.sort_values(inplace=True)
    converted: pd.Dataframe = mg.querymany(genes.values,
                                           scopes='symbol,alias',
                                           species='human',
                                           fields='ensembl.gene,symbol',
                                           as_dataframe=True)
    converted.dropna(axis=0, subset='symbol', inplace=True)
    converted = converted['symbol']
    converted.drop_duplicates(inplace=True)    # Remove duplicated symbols of human genes
    converted = converted[~converted.index.duplicated(keep='first')]    # Remove duplicated symbols of mouse genes
    converted.sort_index(axis=0, inplace=True)
    
    # Remove non converted genes
    for k, v in df_separated.items():
        mask = v.loc[:, f'{k}.names'].isin(converted.index)
        df_separated[k] = v.loc[mask]
        df_separated[k].reset_index(drop=True, inplace=True)
    
    # Fill DEG table
    df_deg: dict[str, list] = {'cell_type': [],
                               'gene': [],
                               'logFC': [],
                               'P.Value': [],
                               'adj.P.Val': []
                               }
    for k, v in df_separated.items():
        df_deg['cell_type'].extend([k] * v.shape[0])
        df_deg['gene'].extend(v.loc[:, f'{k}.names'])
        df_deg['logFC'].extend(v.loc[:, f'{k}.logfoldchanges'])
        df_deg['P.Value'].extend(v.loc[:, f'{k}.pvals'])
        df_deg['adj.P.Val'].extend(v.loc[:, f'{k}.pvals_adj'])
    df_deg = pd.DataFrame(df_deg)
    
    print("Save deg file...")
    deg_file_path = f"{cellphonedb}/{out_file}"
    print(deg_file_path)
    df_deg.to_csv(deg_file_path, sep='\t', index=False)
    
    return deg_file_path


def start() -> None:
    import gc
    import os
    import mygene
    from cellphonedb.utils import db_utils, db_releases_utils
    from cellphonedb.src.core.methods import cpdb_analysis_method, cpdb_statistical_analysis_method, cpdb_degs_analysis_method
    
    # Download database if not existing
    db_ver = "v4.1.0" #4.1.0
    dest = f"{cpdb_dir}/{db_ver}"
    if not os.path.exists(dest):
        db_utils.download_database(dest, db_ver)

    # Load merged data
    dest = f"{checkpoint_dir}/adata_final_merged_raw_norm_annot_nona.h5ad"
    if os.path.exists(dest):
        print("Load merged data...")
        print(dest)
        adata = sc.read_h5ad(dest)
        print(adata)
    else:
        return
    
    # Convert mouse genes to human
    print("\nConvert mouse genes to human")
    mg = mygene.MyGeneInfo()
    converted: pd.Dataframe = mg.querymany(adata.var_names, scopes='symbol,alias', species='human', fields='ensembl.gene,symbol', as_dataframe=True)
    converted.dropna(axis=0, subset='symbol', inplace=True)
    converted = converted['symbol']
    converted.drop_duplicates(inplace=True)    # Remove duplicated symbols of human genes
    converted = converted[~converted.index.duplicated(keep='first')]    # Remove duplicated symbols of mouse genes
    
    # Rename genes
    print("\nRename genes...")
    adata2 = adata.copy()
    adata2 = adata2[:, converted.index]
    adata2.var_names = converted.loc[adata2.var_names]
    
    # Save renamed adata
    print("Save renamed AnnData...")
    dest = f"{checkpoint_dir}/adata_final_merged_raw_norm_hs_names_nona.h5ad"
    print(dest)
    adata2.write_h5ad(dest, compression='gzip')
    
    # CellPhoneDB (cpdb_analysis_method)
    print("\nCellPhoneDB...")
    # Define paths
    #cpdb_file_path: (mandatory) path to the database cellphonedb.zip.
    #meta_file_path: (mandatory) path to the meta file linking cell barcodes to cluster labels metadata.tsv.
    #counts_file_path: (mandatory) paths to normalized counts file (not z-transformed), either in text format or h5ad (recommended) normalised_log_counts.h5ad.
    cpdb_file_path = f"{cpdb_dir}/{db_ver}/cellphonedb.zip"
    counts_file_path = f"{checkpoint_dir}/adata_final_merged_raw_norm_hs_names_nona.h5ad"
    out_path = cellphonedb
    
    # Build meta files
    print("Build meta files...")
    metafile_all = build_meta_file(adata2, cell_type="leiden_merge")   # All cells

    # Free memory
    del adata
    del adata2
    gc.collect()
    
    ### Simple analysis
    # This method will calculate the mean expression of the interacting partners (proteins participating in the interaction)
    # that are expressed in more than threshold percent of cells at each cluster.
    # CellPhoneDB (cpdb_analysis_method)
    if simple_analysis:
        means, deconvoluted = cpdb_analysis_method.call(
            cpdb_file_path=cpdb_file_path,           # mandatory: CellPhoneDB database zip file.
            meta_file_path=metafile_all,             # mandatory: tsv file defining barcodes to cell label.
            counts_file_path=counts_file_path,       # mandatory: normalized count matrix.
            counts_data='hgnc_symbol',               # defines the gene annotation in counts matrix.
            output_path=out_path,                    # Path to save results
            separator='|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
            threshold=0.2,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
            result_precision=4,                      # Sets the rounding for the mean values in significan_means.
            debug=False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
            output_suffix="final_merged_nona"             # Replaces the timestamp in the output files by a user defined string in the  (default: None)
        )
    
    ### Statistical analysis
    # This method will retrieve interactions where the mean expression of the interacting partners
    # (proteins participating in the interaction) displays significant cell state specificity by
    # employing a random shuffling methodology.
    if statistical_analysis:
        # All cells
        deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
            cpdb_file_path=cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
            meta_file_path=metafile_all,                   # mandatory: tsv file defining barcodes to cell label.
            counts_file_path=counts_file_path,             # mandatory: normalized count matrix.
            counts_data='hgnc_symbol',                     # defines the gene annotation in counts matrix.
            iterations=1000,                               # denotes the number of shufflings performed in the analysis.
            threshold=0.2,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
            threads=os.cpu_count() - 1,                    # number of threads to use in the analysis.
            debug_seed=42,                                 # debug randome seed. To disable >=0.
            result_precision=4,                            # Sets the rounding for the mean values in significan_means.
            pvalue=0.05,                                   # P-value threshold to employ for significance.
            subsampling=False,                             # To enable subsampling the data (geometri sketching).
            subsampling_log=False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
            subsampling_num_pc=100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
            subsampling_num_cells=None,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
            separator='|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
            debug=False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
            output_path=out_path,                          # Path to save results.
            output_suffix="final_merged_nona"                   # Replaces the timestamp in the output files by a user defined string in the  (default: None).
        )
    
    ### Differential expression analysis
    # This method will retrieve interactions where at least one of the interacting partners
    # (genes involved in the interaction) is differentially expressed.
    if deg_analysis:
        # For All cells
        degfile_all = build_deg_file(f"{marker_genes_dir}/adata_final_merged_marker_genes_leiden_merge_nona.csv",
                                     "adata_final_merged_deg_nona.tsv")
        deconvoluted, means, relevant_interactions, significant_means = cpdb_degs_analysis_method.call(
            cpdb_file_path=cpdb_file_path,                            # mandatory: CellPhoneDB database zip file.
            meta_file_path=metafile_all,                              # mandatory: tsv file defining barcodes to cell label.
            counts_file_path=counts_file_path,                        # mandatory: normalized count matrix.
            degs_file_path=degfile_all,                               # mandatory: tsv file with DEG to account.
            counts_data='hgnc_symbol',                                # defines the gene annotation in counts matrix.
            threshold=0.2,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
            result_precision=4,                                       # Sets the rounding for the mean values in significan_means.
            separator='|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
            debug=False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
            output_path=cellphonedb,                                  # Path to save results
            output_suffix="final_merged_nona"                              # Replaces the timestamp in the output files by a user defined string in the  (default: None)
        )


start()

print("\n********\n* DONE *\n********")
