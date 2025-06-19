# Perform CellPhoneDB analysis - LeonorSaude10x
#
# Follwing these tutorials
#
# Daniel Ribeiro, Gonçalo Alves 2025
import gc
import os

import pandas as pd

statistical_analysis = True

cpdb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/database"
cellphonedb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb"
cellphonedb_dir_out = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary/summary2"

control = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_significant_means_final_merged_uninjured_nona.txt"
injured_15 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_significant_means_final_merged_injured_15_nona.txt"
injured_60 = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/statistical_analysis_significant_means_final_merged_injured_60_nona.txt"

lineage_colors = {
    'Neuron': 'darkorchid',
    'Meningeal_Vascular': 'slategrey',
    'Immune': 'lime'
}

datset_names = {'Neu': 'Neuron',
                'MeV': 'Meningeal_Vascular',
                'Imm': 'Immune'
                }

def get_cell_types(cci: str,
                   major: bool      # Only major cell types, no clusters
                   
                   ) -> tuple[str, ...]:
    cell_types = cci.split('|')
    if major:
        cell_types[0] = cell_types[0].split('.')[0]
        cell_types[1] = cell_types[1].split('.')[0]

    
    return tuple(cell_types)

def simple_uniprot_to_gene_name(df_simple: pd.DataFrame,      # gene_input.csv dataframe
                                uniprot: str) -> str:
    mask = df_simple.loc[:, "uniprot"] == uniprot
    subset = df_simple[mask]
    # Get fisrt row
    row = subset.iloc[0, :]
    
    return row["gene_name"]


def complex_uniprot_to_gene_name(df_complex: pd.DataFrame,      # complex_input.csv dataframe
                                 df_simple: pd.DataFrame,       # gene_input.csv dataframe
                                 complex_name: str) -> list[str]:
    result = []
    mask = df_complex.loc[:, "complex_name"] == complex_name
    subset = df_complex[mask]
    row = subset.iloc[0, :]
    # Get all possilbe Uniprot ids
    row = row.loc[["uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4", "uniprot_5"]]
    uniprot = []
    uniprot.extend(row[pd.notna(row)].values.tolist())
    for uni in uniprot:
        result.append(simple_uniprot_to_gene_name(df_simple, uni))

    return result


def genes_from_partner(df_complex: pd.DataFrame,      # complex_input.csv dataframe
                       df_simple: pd.DataFrame,       # gene_input.csv dataframe
                       partner: str):
    result = []
    if partner.find("simple:") != -1:
        result.append(simple_uniprot_to_gene_name(df_simple, partner.split(':')[1]))
    elif partner.find("complex:") != -1:
        result.extend(complex_uniprot_to_gene_name(df_complex, df_simple, partner.split(':')[1]))
    
    
    return result


def collect_partners_mp(df_stat: pd.DataFrame,      # statistical_analysis_significant_means dataframe
                        df_complex: pd.DataFrame,   # complex_input.csv dataframe
                        df_simple: pd.DataFrame,    # gene_input.csv dataframe
                        interval: range,
                        unique: bool,               # return only unique names
                        ) -> dict:
    # Progress
    start = min(interval) - 2
    stop = max(interval) - 2
    print(f"Genes for [{start}-{stop}] ccis")
    
    # Columns with CCI
    cols_cci = df_stat.columns[df_stat.columns.str.contains('|', regex=False)]
    df_partner = df_stat.loc[:, ["partner_a", "partner_b"] + cols_cci.to_list()]
    cluster_cci = {}
    for j in interval:
        mask = pd.notna(df_partner.iloc[:, j])
        partner_a = df_partner.loc[mask].loc[:, "partner_a"].to_list()
        partner_b = df_partner.loc[mask].loc[:, "partner_b"].to_list()
        # Solve partner_a genes
        genes_a = []
        for p in partner_a:
            genes_a.extend(genes_from_partner(df_complex, df_simple, p))
        genes_a.sort()
        # Solve partner_b genes
        genes_b = []
        for p in partner_b:
            genes_b.extend(genes_from_partner(df_complex, df_simple, p))
        genes_b.sort()
        # Join partner genes
        genes_joint = genes_a + genes_b
        if unique:
            genes_joint = list(set(genes_joint))    # Unique names
        genes_joint.sort()
        cluster_cci[df_partner.columns[j]] = genes_joint
            
    gc.collect()
    return cluster_cci


def collect_partners(df_stat: pd.DataFrame,      # statistical_analysis_significant_means dataframe
                     df_complex: pd.DataFrame,   # complex_input.csv dataframe
                     df_simple: pd.DataFrame,    # gene_input.csv dataframe
                     unique: bool = False,       # return only unique gene names per partner group
                     n_proc: int = None) -> dict:
    from multiprocessing.pool import Pool
    
    # Columns with CCI
    cols_cci = df_stat.columns[df_stat.columns.str.contains('|', regex=False)]
    df_partner = df_stat.loc[:, ["partner_a", "partner_b"] + cols_cci.to_list()]
    steps = 1000  # Count in steps of 1000
    # Create list of ranges
    start = 2
    stop = len(df_partner.columns)
    end = start + steps
    ranges = []
    while start < stop:
        ranges.append(range(start, end))
        start += steps
        end = end + steps if (end + steps) <= stop else stop
    
    cluster_cci = {}
    with Pool(processes=n_proc) as p:
        # Prepare args for batches of 1000
        args = [(df_stat, df_complex, df_simple, interval, unique)
                for interval in ranges]
        # Prepare args for within comparisons
        dicts = list(p.starmap(collect_partners_mp, args))
        for d in dicts:
            if d is None:
                continue
            cluster_cci.update(d)
        
        # Free memory
        gc.collect()

    #import sys
    #    # Progress
    #    if (counter % steps == 0) or (counter == (len(df_partner.columns) - 1)):
    #        if counter != 0:
    #            sys.stdout.write("\r\033[2K")
    #        print(f"Genes for: {counter}/{total_cci} ccis", end="", flush=True)
    #    counter += 1
    #
    #print()
    

    return cluster_cci


def get_source_targets(summary_df: pd.DataFrame,
                       major_cell_types: bool,      # Only use major cell types, no clusters
                       unique_cci: bool,            # Only report unique cci. Useful when collapsing cell types
                       ) -> tuple[pd.DataFrame, pd.DataFrame]:
    from itertools import product
    import numpy as np
    # Build a source/target matrix - row = source, col = target
    # Get all possible sources and targets
    sources = []
    targets = []
    for row in summary_df.index:
        cci = get_cell_types(row, major_cell_types)
        sources.append(cci[0])
        targets.append(cci[1])
    sources = list(set(sources))
    sources.sort()
    targets = list(set(targets))
    targets.sort()
    mtx_source_target = pd.DataFrame(np.zeros((len(sources), len(targets))), index=sources, columns=targets, dtype='int64')
    
    # For every interaction add a +1 weight
    # Get all ccis for each src/tgt combination
    keys = list(product(mtx_source_target.index, mtx_source_target.columns))
    interactions: dict[str, list] = {k: [] for k in keys}
    for row in summary_df.index:
        ccis = summary_df.loc[row, summary_df.loc[row, :].notna()].values.tolist()
        src_tgt = get_cell_types(row, major_cell_types)
        interactions[src_tgt].extend(ccis)
    
    # Only collect unique cci
    if unique_cci:
        for row in summary_df.index:
            src_tgt = get_cell_types(row, major_cell_types)
            interactions[src_tgt] = list(set(interactions[src_tgt]))
            interactions[src_tgt].sort()
            
    # Add weights
    for src_tgt in interactions:
        mtx_source_target.loc[src_tgt[0], src_tgt[1]] = len(interactions[src_tgt])
    
    # holoviews requires input data containing the following columns:
    # “source”, “target”, “weight”
    df_source_target = {"source": [],
                        "target": [],
                        "value": []}
     
    # Get indexes for src_tgt
    nodes = [src_tgt[0] if not major_cell_types else src_tgt[0].split('.')[0] for src_tgt in interactions]
    nodes.extend([src_tgt[1] if not major_cell_types else src_tgt[1].split('.')[0] for src_tgt in interactions])
    nodes = list(set(nodes))
    nodes.sort()
    nodes = {n: i for i, n in enumerate(nodes)}
    # print("nodes.....")
    # print(nodes)
    for src_tgt in interactions:
        src = src_tgt[0]
        tgt = src_tgt[1]
        df_source_target["source"].append(nodes[src])
        df_source_target["target"].append(nodes[tgt])
        df_source_target["value"].append(mtx_source_target.loc[src, tgt])
    
    df_source_target = pd.DataFrame(df_source_target)
    # nodes = {"index": nodes.values(),
    #          "names": nodes.keys(),
    #          "group": [0] * len(nodes.values())}
    # nodes = pd.DataFrame(nodes)

    # Construct nodes correctly
    node_list = list(nodes.keys())  # Extract node names
    node_indices = list(nodes.values())  # Extract corresponding indices
    
    nodes = pd.DataFrame({
        "index": node_indices,
        "names": node_list,
        "group": [0] * len(node_list)   
    })

    
    # Also return index names

    print("hey")
    print(df_source_target)
    print("2ndhey")
    print(nodes)    

    return df_source_target, nodes
    





def start(n_proc: int = None) -> None:
    import pandas as pd

    ### Define multiple input files
    input_files = {
        "control": control,
        "injured_15": injured_15,
        "injured_60": injured_60
    }

    ### Statistical analysis - For each condition
    for condition, filepath in input_files.items():
        if not os.path.exists(filepath):
            print(f"File not found for condition '{condition}': {filepath}")
            continue

        print(f"\nProcessing condition: {condition}")
        df = pd.read_csv(filepath, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
        
        # Columns with CCI
        cols_cci = df.columns[df.columns.str.contains('|', regex=False)]
        df_cci = df.loc[:, ["id_cp_interaction", "interacting_pair"] + cols_cci.to_list()]
        new_index = []
        for i in df_cci.index:
            new_index.append(f"{df_cci.loc[i, 'interacting_pair']}_{df_cci.loc[i, 'id_cp_interaction']}")
        df_cci.index = new_index
        if df_cci.index.has_duplicates:
            raise ValueError(f"Index has duplicates in {condition}!")
        df_cci.index.name = 'cci'
        df_cci.drop(labels=["id_cp_interaction", "interacting_pair"], axis=1, inplace=True)

        # Summarize CCIs
        cluster_cci = {}
        for col in df_cci.columns:
            cluster_cci[col] = []
        for j in range(len(df_cci.columns)):
            mask = pd.isna(df_cci.iloc[:, j])
            cluster_cci[df_cci.columns[j]].extend(df_cci.iloc[:, j].loc[~mask].index.to_list())
        for col in df_cci.columns:
            cluster_cci[col].sort()

        # Save summary of CCIs
        dest = f"{cellphonedb_dir_out}/summary_significant_cci_means_{condition}.txt"
        cluster_cci = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cluster_cci.items()]))
        cluster_cci.T.to_csv(dest, index=True, header=False, sep='\t')
        
        # Collect genes for interactions
        print(f"Finding genes for {condition}")
        gene_input = pd.read_csv(f"{cpdb_dir}/v4.1.0/gene_input.csv", sep=',', header=0)
        complex_input = pd.read_csv(f"{cpdb_dir}/v4.1.0/complex_input.csv", sep=',', header=0)
        cci_genes = collect_partners(df_stat=df, df_complex=complex_input, df_simple=gene_input, n_proc=n_proc)
        
        # Save genes
        dest = f"{cellphonedb_dir_out}/summary_significant_cci_genes_{condition}.txt"
        cci_genes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cci_genes.items()]))
        cci_genes.T.to_csv(dest, index=True, header=False, sep='\t')




# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=os.cpu_count() - 1)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
    