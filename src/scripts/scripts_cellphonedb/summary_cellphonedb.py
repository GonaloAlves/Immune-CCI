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
cellphonedb_dir = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/excels/edge_list"
cellphonedb_dir_out = "/home/makowlg/Documents/Immune-CCI/src/cellphonedb/summary"


datset_names = {'Neu': 'Neuron',
                'MeV': 'Meningeal_Vascular',
                'Imm': 'Immune'
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


    return cluster_cci

def get_genes_from_id_cp_interaction(interaction_id: str, df_stat: pd.DataFrame, df_complex: pd.DataFrame, df_simple: pd.DataFrame) -> list[str]:
    
    # Find row in df_stat for this interaction
    row = df_stat[df_stat['id_cp_interaction'] == interaction_id]
    if row.empty:
        return []

    partner_a = row['partner_a'].values[0]
    partner_b = row['partner_b'].values[0]
    
    genes = []
    # Handle partner_a
    genes += genes_from_partner(df_complex, df_simple, partner_a)
    # Handle partner_b
    genes += genes_from_partner(df_complex, df_simple, partner_b)

    return sorted(set(genes))  # optional: remove duplicates


def enrich_simplified_excel_with_genes( simplified_excel_path: str, output_excel_path: str, cpdb_gene_input_path: str, cpdb_complex_input_path: str, cpdb_stat_path: str):

    # Load CellPhoneDB inputs
    df_simple = pd.read_csv(cpdb_gene_input_path, sep=',')
    df_complex = pd.read_csv(cpdb_complex_input_path, sep=',')
    df_stat = pd.read_csv(cpdb_stat_path, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})

    # Load all sheets
    xl = pd.read_excel(simplified_excel_path, sheet_name=None)

    updated_sheets = {}
    for sheet_name, df in xl.items():
        print(f"Processing {sheet_name}...")
        genes_list = []

        for interaction_id in df["id_cp_interaction"]:
            genes = get_genes_from_id_cp_interaction(interaction_id, df_stat, df_complex, df_simple)
            genes_list.append(", ".join(genes))  # Join as string for Excel

        df["genes"] = genes_list
        updated_sheets[sheet_name] = df

    # Write to new Excel
    with pd.ExcelWriter(output_excel_path, engine='xlsxwriter') as writer:
        for sheet_name, df in updated_sheets.items():
            df.to_excel(writer, index=False, sheet_name=sheet_name)

    print("✅ Done writing enriched Excel.")

def start(n_proc: int = None) -> None:
    import pandas as pd

    # ### Statistical analysis - Summarize by cci
    # if statistical_analysis:
    #     # Load the multiple statistical analysis
    #     dest = f"{cellphonedb_dir}/statistical_analysis_significant_means_final_merged.txt"
    #     print("Summarizing ", dest)
    #     df = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"})
    #     # Columns with CCI
    #     cols_cci = df.columns[df.columns.str.contains('|', regex=False)]
    #     df_cci = df.loc[:, ["id_cp_interaction", "interacting_pair"] + cols_cci.to_list()]
    #     # There are duplicate names in interacting_pair, append id to interaction
    #     new_index = []
    #     for i in df_cci.index:
    #         new_index.append(f"{df_cci.loc[i, 'interacting_pair']}_{df_cci.loc[i, 'id_cp_interaction']}")
    #     df_cci.index = new_index
    #     if df_cci.index.has_duplicates:
    #         raise ValueError("Index has duplicates!")
    #     df_cci.index.name = 'cci'
    #     df_cci.drop(labels=["id_cp_interaction", "interacting_pair"], axis=1, inplace=True)

    #     # Collect cci for each interaction
    #     cluster_cci = {}
    #     for col in df_cci.columns:
    #         cluster_cci[col] = []
    #     for j in range(len(df_cci.columns)):
    #         mask = pd.isna(df_cci.iloc[:, j])
    #         cluster_cci[df_cci.columns[j]].extend(df_cci.iloc[:, j].loc[~mask].index.to_list())
    #     # Sort cci names
    #     for col in df_cci.columns:
    #         cluster_cci[col].sort()

    #     # Save
    #     dest = f"{cellphonedb_dir_out}/summary_significant_cci_means_final_merged.txt"
    #     cluster_cci = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cluster_cci.items()]))
    #     cluster_cci.T.to_csv(dest, index=True, header=False, sep='\t')
    
    ### Statistical analysis - Summarize by cci genes
    # Collect genes from each partner pair
    if statistical_analysis:
        gene_input = pd.read_csv(f"{cpdb_dir}/v4.1.0/gene_input.csv", sep=',', header=0)
        complex_input = pd.read_csv(f"{cpdb_dir}/v4.1.0/complex_input.csv", sep=',', header=0)
        # Find genes for each interaction
        
        dest = f"{cellphonedb_dir}/statistical_analysis_significant_means_final_merged.txt"
        print("\nFinding genes for", dest)
        df_stat = pd.read_csv(dest, sep='\t', dtype={"gene_a": "string", "gene_b": "string"}, low_memory=False)
        cci_genes = collect_partners(df_stat=df_stat, df_complex=complex_input, df_simple=gene_input)
        # Save
        dest = f"{cellphonedb_dir_out}/summary_significant_cci_genes_final_merged.txt"
        cci_genes = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cci_genes.items()]))
        cci_genes.T.to_csv(dest, index=True, header=False, sep='\t')
##
        enrich_simplified_excel_with_genes(
        simplified_excel_path="/path/to/simplified_significant_interactions_15.xlsx",
        output_excel_path="/path/to/enriched_significant_interactions_15_with_genes.xlsx",
        cpdb_gene_input_path="/path/to/gene_input.csv",
        cpdb_complex_input_path="/path/to/complex_input.csv",
        cpdb_stat_path="/path/to/statistical_analysis_significant_means_final_merged.txt"
        )


    

# main guard required because processes are spawn (compatible with Windows)
if __name__ == '__main__':

    import multiprocessing as mp
    try:
        mp.set_start_method('spawn', force=True)   # Ensure Windows compatibility
        start(n_proc=os.cpu_count() - 1)

        print("\n********\n* DONE *\n********")
    except RuntimeError:
        raise
    