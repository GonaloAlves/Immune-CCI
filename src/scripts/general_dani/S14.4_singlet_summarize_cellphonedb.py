# Perform CellPhoneDB analysis - LeonorSaude10x
#
# Follwing these tutorials
# https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T01_Method1.ipynb
# https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T01_Method2.ipynb
#
# Daniel Ribeiro, 2023
import gc
import os

import pandas as pd

statistical_analysis = True
datset_names = {'Neu': 'Neuron',
                'Oli': 'Oligodendrocyte',
                'Ast': 'Astrocyte',
                'MeV': 'Meningeal_Vascular',
                'Imm': 'Immune',
                }


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
    from src.helper_functions import get_cell_types
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
    for src_tgt in interactions:
        src = src_tgt[0]
        tgt = src_tgt[1]
        df_source_target["source"].append(nodes[src])
        df_source_target["target"].append(nodes[tgt])
        df_source_target["value"].append(mtx_source_target.loc[src, tgt])
    
    df_source_target = pd.DataFrame(df_source_target)
    nodes = {"index": nodes.values(),
             "names": nodes.keys(),
             "group": [0] * len(nodes.values())}
    nodes = pd.DataFrame(nodes)
    
    # Also return index names
    return df_source_target, nodes
    

def chord_diagram(chord_input: pd.DataFrame,
                  node_info: pd.DataFrame,
                  out_file: str) -> None:
    import holoviews as hv
    from bokeh.io import export_svgs, export_png
    from holoviews import opts, dim
    from matplotlib.colors import ListedColormap
    from src.globals import lineage_colors
    
    # Initialize
    hv.extension('bokeh', 'matplotlib')       # Bokeh backend to draw chords
    hv.output(size=500, dpi=300)
    
    # Get colors
    cmap = list(set(node_info["names"].values.tolist()))
    cmap.sort()
    # Ensure that lineages are ok
    for i in range(len(cmap)):
        if '.' in cmap[i]:
            for k in datset_names:
                if cmap[i].startswith(k):
                    cmap[i] = datset_names[k]
                    break
        else:
            cmap[i] = datset_names[cmap[i]]
    cmap = ListedColormap([lineage_colors[c] for c in cmap], len(cmap))
    # Build chord diagram
    node_info = hv.Dataset(node_info, "index")
    chord = hv.Chord(data=(chord_input, node_info))
    plot = chord.opts(opts.Chord(cmap=cmap, edge_cmap=cmap, edge_color=dim('source').str(),
                      labels='names', node_color=dim('index').str(), label_text_font_size='24pt'))
    
    p = hv.render(plot, backend='bokeh')
    #p.output_backend = "svg"
    #export_svgs(p, filename=out_file)
    export_png(plot, filename=out_file)


def gprofiler_analysis(df_genes: pd.DataFrame,
                       major_cell_types: bool,
                       unique_genes: bool) -> dict:
    from gprofiler import GProfiler
    from itertools import product
    from src.helper_functions import get_cell_types
    
    # Get cell types
    print("Collecting genes...")
    sources = []
    targets = []
    for row in df_genes.index:
        cci = get_cell_types(row, major_cell_types)
        sources.append(cci[0])
        targets.append(cci[1])
    sources = list(set(sources))
    sources.sort()
    targets = list(set(targets))
    targets.sort()
    keys = list(product(sources, targets))
    interactions = {k: [] for k in keys}
    
    # Make faster preventing preventing memory realloc
    sizes = {k: 0 for k in keys}
    n = df_genes.notna().sum(axis=1)
    for row in df_genes.index:
        src_tgt = get_cell_types(row, major_cell_types)
        sizes[src_tgt] += n.loc[row]
    for k in sizes:
        interactions[k] = [None] * sizes[k]
    
    # Fill interactions
    pos = {k: 0 for k in keys}
    for row in df_genes.index:
        src_tgt = get_cell_types(row, major_cell_types)
        genes = df_genes.loc[row, df_genes.loc[row, :].notna()].values.tolist()
        i = 0
        while (i < len(genes)) and (pos[src_tgt] < len(interactions[src_tgt])):
            interactions[src_tgt][pos[src_tgt] + i] = genes[i]
            i += 1
        pos[src_tgt] += len(genes)
    
    # Only collect unique genes
    if unique_genes:
        for row in df_genes.index:
            src_tgt = get_cell_types(row, major_cell_types)
            interactions[src_tgt] = list(set(interactions[src_tgt]))
            interactions[src_tgt].sort()
    interactions = pd.DataFrame(dict([('|'.join(k), pd.Series(v)) for k, v in interactions.items()])).T
    
    # Initialize and run gProfiler analysis
    gp = GProfiler(return_dataframe=True)
    # Run gProfiler
    print("Running gProfiler...")
    queries = {row: interactions.loc[row, interactions.loc[row, :].notna()].values.tolist() for row in interactions.index}
    results = {}
    for q in queries:
        results[q] = gp.profile(organism='hsapiens',    # Because genes are converted to hsapiens
                                query=queries[q],
                                sources=['GO:BP', 'KEGG', 'REAC', 'TF'],
                                no_iea=True)
        results[q].sort_values(by=['source', 'p_value'], ascending=[True, True], inplace=True)
        results[q].drop(columns='description', axis=1, inplace=True)
    
    return results


def start(n_proc: int = None) -> None:
    import pandas as pd
    import src.globals   # noqa:F401
    from src.globals import cellphonedb_dir, cpdb_dir, gprofiler_dir

    #injury_day = ["uninjured.0", "sham.15", "injured.15", "injured.60"]
    injury_day = ["sham.15", "injured.15", "injured.60"]
    ### Statistical analysis - Summarize by cci
    if statistical_analysis:
        # Load the multiple statistical analysis
        for inj in injury_day:
            inj2 = inj.replace('.', '_')
            dest = f"{cellphonedb_dir}/statistical_analysis_significant_means_final_merged_{inj2}.txt"
            print("Summarizing ", dest)
            df = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
            # Columns with CCI
            cols_cci = df.columns[df.columns.str.contains('|', regex=False)]
            df_cci = df.loc[:, ["id_cp_interaction", "interacting_pair"] + cols_cci.to_list()]
            # There are duplicate names in interacting_pair, append id to interaction
            new_index = []
            for i in df_cci.index:
                new_index.append(f"{df_cci.loc[i, 'interacting_pair']}_{df_cci.loc[i, 'id_cp_interaction']}")
            df_cci.index = new_index
            if df_cci.index.has_duplicates:
                raise ValueError("Index has duplicates!")
            df_cci.index.name = 'cci'
            df_cci.drop(labels=["id_cp_interaction", "interacting_pair"], axis=1, inplace=True)
            # Collect cci for each interaction
            cluster_cci = {}
            for col in df_cci.columns:
                cluster_cci[col] = []
            for j in range(len(df_cci.columns)):
                mask = pd.isna(df_cci.iloc[:, j])
                cluster_cci[df_cci.columns[j]].extend(df_cci.iloc[:, j].loc[~mask].index.to_list())
            # Sort cci names
            for col in df_cci.columns:
                cluster_cci[col].sort()
    
            # Save
            dest = f"{cellphonedb_dir}/summary_significant_cci_means_final_merged_{inj2}.txt"
            cluster_cci = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cluster_cci.items()]))
            cluster_cci.T.to_csv(dest, index=True, header=False, sep='\t')
    
    ### Statistical analysis - Summarize by cci genes
    # Collect genes from each partner pair
    if statistical_analysis:
        gene_input = pd.read_csv(f"{cpdb_dir}/v4.1.0/gene_input.csv", sep=',', header=0)
        complex_input = pd.read_csv(f"{cpdb_dir}/v4.1.0/complex_input.csv", sep=',', header=0)
        # Find genes for each interaction
        for inj in injury_day:
            inj2 = inj.replace('.', '_')
            dest = f"{cellphonedb_dir}/statistical_analysis_significant_means_final_merged_{inj2}.txt"
            print("\nFinding genes for", dest)
            df_stat = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
            cci_genes = collect_partners(df_stat=df_stat, df_complex=complex_input, df_simple=gene_input)
            # Save
            dest = f"{cellphonedb_dir}/summary_significant_cci_genes_final_merged_{inj2}.txt"
            cci_genes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cci_genes.items()]))
            cci_genes.T.to_csv(dest, index=True, header=False, sep='\t')

    ## Build chord diagram for cci and cci genes
    if statistical_analysis:
        # Chord for major cell types
        for inj in injury_day:
            inj2 = inj.replace('.', '_')
            dest = f"{cellphonedb_dir}/summary_significant_means_final_merged_{inj2}.txt"
            print(f"\nBuild chord diagram for {dest}")
            df = pd.read_csv(dest, index_col=0, header=None, sep='\t', low_memory=False)
            
            # Unique interactions
            chord_input, node_info = get_source_targets(df, major_cell_types=True, unique_cci=True)
            dest = f"{cellphonedb_dir}/chord_significant_final_merged_unique_cci_{inj2}.png"
            chord_diagram(chord_input, node_info, dest)
            print(dest)
            # Save matrices
            dest = f"{cellphonedb_dir}/chord_input_major_cell_types_unique_cci_{inj2}.txt"
            chord_input.to_csv(dest, sep='\t')
            print(dest)
            dest = f"{cellphonedb_dir}/node_info_major_cell_types_unique_cci_{inj2}.txt"
            node_info.to_csv(dest, sep='\t')
            print(dest)
                      
            # All interactions
            chord_input, node_info = get_source_targets(df, major_cell_types=True, unique_cci=False)
            dest = f"{cellphonedb_dir}/chord_significant_final_merged_all_cci_{inj2}.png"
            chord_diagram(chord_input, node_info, dest)
            print(dest)
            # Save matrices
            dest = f"{cellphonedb_dir}/chord_input_major_cell_types_all_cci_{inj2}.txt"
            chord_input.to_csv(dest, sep='\t')
            print(dest)
            dest = f"{cellphonedb_dir}/node_info_major_cell_types_all_cci_{inj2}.txt"
            node_info.to_csv(dest, sep='\t')
            print(dest)
    
        # Chord for clusters
        for inj in injury_day:
            inj2 = inj.replace('.', '_')
            dest = f"{cellphonedb_dir}/summary_significant_means_final_merged_{inj2}.txt"
            print(f"\nBuild chord diagram for {dest}")
            df = pd.read_csv(dest, index_col=0, header=None, sep='\t', low_memory=False)
            chord_input, node_info = get_source_targets(df, major_cell_types=False, unique_cci=True)
            dest = f"{cellphonedb_dir}/chord_significant_final_merged_unique_cci_clusters_{inj2}.png"
            chord_diagram(chord_input, node_info, out_file=dest)
    
    ## gProfiler analysis
    # Do gProfiler analysis of terms
    if statistical_analysis:
        for inj in injury_day:
            # Find genes for each interaction
            inj2 = inj.replace('.', '_')
            dest = f"{cellphonedb_dir}/summary_significant_cci_genes_final_merged_{inj2}.txt"
            print("\ngProfiler analysis for", dest)
            df_genes = pd.read_csv(dest, sep='\t', index_col=0, header=None, dtype='string', low_memory=False)
            gprofiler_res = gprofiler_analysis(df_genes=df_genes, major_cell_types=True, unique_genes=True)
            
            # Save
            for r in gprofiler_res:
                dest = f"{gprofiler_dir}/gProfiler.gGOSt_cci_unique_genes_final_merged_{r}_{inj2}.txt"
                print(f"Save {dest}")
                gprofiler_res[r].to_csv(dest, index=True, header=True, sep='\t')
    

# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=os.cpu_count() - 1)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
