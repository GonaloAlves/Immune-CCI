# # Analysis of all samples - LeonorSaude10x
# The following steps are typically performed to analyze a scRNA-seq sample:
#
# Global variables for the analysis.
# Setup folders and modules.
#
#   Daniel Ribeiro, 2022
#
import scanpy as sc
import warnings

warnings.simplefilter("ignore", UserWarning)

# Dataset suffix
f_suffix = 'automax'
cell_suffix = '0k'

# Downsampling
downsample = False
downsample_target = 3500

# Quality control variables
min_counts = 430
max_counts = 6000
min_genes = 480
max_genes = 4000
percent_mito = 1
min_cells = 5
doublet_removal = False

# Batch correction
# NOTE: Some batch algorithms are too strong and may distort the data.
#batch_combat = False
batch_mnn = True
batch_scanorama = True
load_batch = True

# High variable genes
n_hvg = 2000
n_hvg_clean = 4000
n_hvg_spi = 4000
hvg_flavor = 'seurat'

# Clustering
# Try 3 n_neighbours = 5, 15 (default), 30
#n_neighbors = [5, 15, 30]
n_neighbors = [30]
n_neighbors_final = [15]
resolution = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]
resolution_tests = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Directories
cell_cycle_dir = ""
cell_lineage_dir = ""
cluster_fusion_dir = ""
cluster_identity_dir = ""
cpdb_dir = ""
input_dir = ""
external_dir = ""
gene_sets_dir = ""
output_dir = ""
checkpoint_dir = ""
preprocessing_dir = ""
normalization_dir = ""
feature_dir = ""
dim_reduction_dir = ""
clustering_dir = ""
confounding_dir = ""
lineage_doublets_dir = ""
lineage_marker_dir = ""
putative_lineage_dir = ""
dotplot_dir = ""
scatterplot_dir = ""
violinplot_dir = ""
hca_dir = ""
barplot_dir = ""
venn_diagram_dir = ""
marker_genes_dir = ""
marker_vis_dir = ""
unique_markers_dir = ""
clustree_dir = ""
annotation_dir = ""
scMCA_annotation_dir = ""
fractions_dir = ""
gsea_dir = ""
scdiffcom_dir = ""
cellphonedb_dir = ""
gprofiler_dir = ""
revigo_dir = ""
sc_clustering_dir = ""

# Sample groups
injury_day = ['uninjured_0', 'sham_15', 'injured_15', 'injured_60']
injury_condition = ['uninjured_0_central', 'sham_15_rostral', 'sham_15_caudal', 'injured_15_rostral', 'injured_15_caudal', 'injured_60_rostral', 'injured_60_caudal']
injury_region = ['uninjured_central', 'sham_rostral', 'sham_caudal', 'injured_rostral', 'injured_caudal']
collection_region = ['central', 'rostral', 'caudal']
injury_region_no_central = ['sham_rostral', 'sham_caudal', 'injured_rostral', 'injured_caudal']
collection_region_no_central = ['rostral', 'caudal']

# Datasets
data: dict[str, sc.AnnData] = {}
# Enumerate dataset names and corresponding batch keys for correction
data_sets = {
    'adata': None,
    'adata_seurat_mouse': 'mouse_id',
    'adata_spi': 'mouse_id',
    'adata_final': 'mouse_id',
    'adata_final_merged': 'mouse_id'
}

## Major datasets
datasets_divided = {
    'Neuron': ['Neuron',
               'Neuroepithelial'],
    'Oligodendrocyte': ['Oligodendrocyte', 'OPC', 'Schwann'],
    'Astrocyte': ['Astrocyte'],
    'Meningeal_Vascular': ['Meningeal',
                           'Endothelial',
                           'Fibroblast',
                           'Pericyte',
                           'Stroma',
                           'Muscle',
                           'Glia',          # classified as Müller glia
                           'Epithelial',
                           'Bone',
                           'Erythroid'],
    'Immune': ['Immune'],
}

# Small datasets
datasets_small_names = {
    'Neuron': 'Neu',
    'Oligodendrocyte': 'Oli',
    'Astrocyte': 'Ast',
    'Meningeal_Vascular': 'MeV',
    'Immune': 'Imm'
}
    

# Cell removal
cluster_cell_removal = {
    'adata_seurat_mouse': {0.4: [0, 1]},        # Resolution: Cluster number
    'adata_seurat_mouse_clean': {0.4: 10}       # Min expression: Resolution
}

# Fixed colors
lineage_colors = {
    'Neuron': 'darkorchid',
    'Oligodendrocyte': 'orange',
    'Astrocyte': 'skyblue',
    'Meningeal': 'slategrey',
    'Meningeal_Vascular': 'slategrey',
    'Immune': 'lime',
    'Endothelial': 'crimson',
    'Fibroblast': 'aqua',
    'Pericyte': 'forestgreen',
    'Stroma': 'black',
    'Neuroepithelial': 'mediumblue',
    'Muscle': 'chocolate',
    'Erythroid': 'lightsalmon',     # few cells
    'Bone': 'steelblue',            # few cells
    'Epithelial': 'yellowgreen',    # few cells
    'Glia': 'cornflowerblue'        # few cells, classified as Müller Glia
}

lineage_shades = {
    'Neuron': {'Purples': ['excitatory'],
               'Blues': ['inhibitory'],
               'Greens': ['cholinergic'],
               'Oranges': ['CSF-cN', 'ependymal'],
               'Greys': ['mixed', 'undetermined'],
               },
    'Oligodendrocyte': {},
    'Astrocyte': {},
    'Meningeal_Vascular': {},
    'Immune': {},
}

proportion_colors = {
    'all': '#899499',
    'uninjured_0': '#428bca',
    'sham_15': '#5cb85c',
    'injured_15': '#ff9922',
    'injured_60': '#d9534f',
    'uninjured_0_central': '#428bca',
    'sham_15_rostral': '#5cb85c',
    'sham_15_caudal': '#addbad',
    'injured_15_rostral': '#ff9922',
    'injured_15_caudal': '#f7bf81',
    'injured_60_rostral': '#d9534f',
    'injured_60_caudal': '#eca9a7',
    'caudal': '#b3bfa1',
    'rostral': '#d2a295',
    'central': '#428bca',
    'injured_caudal': '#9b636b',
    'injured_rostral': '#d2a295',
    'sham_caudal': '#5cb85c',
    'sham_rostral': '#addbad',
    'uninjured_central': '#428bca',
}
    

# Lineage resolution tests for clustree
lineage_resolution_tests = {
    'Neuron': [0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    'Oligodendrocyte': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5],
    'Astrocyte': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2],
    'Meningeal_Vascular': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    'Immune': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2],
}
weighted_SC_resolution_tests = {
    'Neuron': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    'Oligodendrocyte': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5],
    'Astrocyte': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2],
    'Meningeal_Vascular': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    'Immune': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2],
}

# Cluster resolutions to analyze in dotplots
lineage_resolution = {
    'Neuron': [5, 10, 13],
    'Oligodendrocyte': [1.4],
    'Astrocyte': [0.7],
    'Meningeal_Vascular': [1.4],
    'Immune': [1]
}
lineage_resolution_uinj = {
    'Neuron_uinj': [0.5, 0.8, 1.5, 5],
    'Oligodendrocyte_uinj': [0.9, 1.4],
    'Astrocyte_uinj': [0.7],
    'Meningeal_Vascular_uinj': [0.1, 1],
    'Immune_uinj': [0.3, 0.4, 0.8]
}
lineage_resolution_final = {
    'Neuron': [5],
    'Oligodendrocyte': [0.9],
    'Astrocyte': [0.7],
    'Meningeal_Vascular': [1.4],
    'Immune': [0.7]
}

# Senescence cluster resolutions
senescence_resolution = {
    'Neuron': {'NES_CELLAGE': [0.7],
               'NES_SENMAYO': [0.4],
               'NES_REACTOME_SASP': [0.5],
               'NES_FRIEDMAN_UP': [0.5],
               'NES_FBR_UP': [0.6],
               'NES_INTEGRATED': [0.4],
               #'GSVA_CELLAGE': [0.4],
               #'GSVA_SENMAYO': [0.4],
               #'GSVA_REACTOME_SASP': [0.5],
               },
    'Oligodendrocyte': {'NES_CELLAGE': [0.8],
                        'NES_SENMAYO': [0.7],
                        'NES_REACTOME_SASP': [0.8],
                        'NES_FRIEDMAN_UP': [0.6],
                        'NES_FBR_UP': [0.5],
                        'NES_INTEGRATED': [0.4],
                        #'GSVA_CELLAGE': [0.9],
                        #'GSVA_SENMAYO': [0.6],
                        #'GSVA_REACTOME_SASP': [0.9],
                        },
    'Astrocyte': {'NES_CELLAGE': [0.7],
                  'NES_SENMAYO': [0.6],
                  'NES_REACTOME_SASP': [0.4],
                  'NES_FRIEDMAN_UP': [0.4],
                  'NES_FBR_UP': [0.6],
                  'NES_INTEGRATED': [0.7],
                  #'GSVA_CELLAGE': [0.6],
                  #'GSVA_SENMAYO': [0.5],
                  #'GSVA_REACTOME_SASP': [0.7],
                  },
    'Meningeal_Vascular': {'NES_CELLAGE': [1.2],
                           'NES_SENMAYO': [1.3],
                           'NES_REACTOME_SASP': [0.8],
                           'NES_FRIEDMAN_UP': [0.9],
                           'NES_FBR_UP': [1.2],
                           'NES_INTEGRATED': [0.8],
                           #'GSVA_CELLAGE': [0.9],
                           #'GSVA_SENMAYO': [0.6],
                           #'GSVA_REACTOME_SASP': [0.3],
                           },
    'Immune': {'NES_CELLAGE': [1.2],
               'NES_SENMAYO': [1.0],
               'NES_REACTOME_SASP': [1.4],
               'NES_FRIEDMAN_UP': [1.0],
               'NES_FBR_UP': [1.3],
               'NES_INTEGRATED': [0.8],
               #'GSVA_CELLAGE': [0.6],
               #'GSVA_SENMAYO': [0.5],
               #'GSVA_REACTOME_SASP': [0.6],
               },
}

# Cluster annotation
cluster_annotation = {
    'Neuron': {'excitatory': 'Excit',
               'inhibitory': 'Inhib',
               'cholinergic': 'Choli',
               'CSF-cN': 'CSF-cN',
               'ependymal': 'Epend',
               'mixed': 'Mix',
               'undetermined': 'NA',
               },
    'Oligodendrocyte': {},
    'Astrocyte': {},
    'Meningeal_Vascular': {},
    'Immune': {},
    

}


# GSEA vars
# Gene sets for GSEApy. name are taken from https://maayanlab.cloud/Enrichr/#libraries
gene_sets = {
    "H:hallmark": ["MSigDB_Hallmark_2020"],
    "C2:curarted": ["KEGG_2016"],
    "C5:ontology": ["GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"],
}
gene_sets_senescence = {
    "FBR_UP": None,
    "FBR_DN": None,
    "FRIEDMAN_UP": None,
    "REACTOME_SASP": None,
    "SENMAYO": None,
    "CELLAGE": None,
    "CELLAGE_UP": None,
    "CELLAGE_DN": None,
}


def init() -> None:
    import os
    
    # Load Scanpy
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.settings.set_figure_params(dpi=300, dpi_save=300, color_map='viridis', facecolor='white')  # Get high quality figs
    sc._settings.ScanpyConfig(n_jobs=os.cpu_count())

    #sc.logging.print_header()

    #print("\nCurrent working directory:")
    print(os.getcwd())

    # Create output dirs
    #global output_dir
    output_dir = "output"
    global cell_cycle_dir
    cell_cycle_dir = "cell_cycle"
    global cell_lineage_dir
    cell_lineage_dir = "cell_lineage"
    global cluster_fusion_dir
    cluster_fusion_dir = "cluster_fusion"
    global cluster_identity_dir
    cluster_identity_dir = "cluster_identity"
    global cpdb_dir
    cpdb_dir = "cpdb"
    global external_dir
    external_dir = "external"
    global gene_sets_dir
    gene_sets_dir = "gene_sets"
    global input_dir
    input_dir = "input"

    #print("Initialize...")
    global checkpoint_dir
    global preprocessing_dir
    global normalization_dir
    global feature_dir
    global dim_reduction_dir
    global clustering_dir
    global confounding_dir
    global lineage_doublets_dir
    global lineage_marker_dir
    global putative_lineage_dir
    global dotplot_dir
    global scatterplot_dir
    global violinplot_dir
    global hca_dir
    global barplot_dir
    global venn_diagram_dir
    global marker_genes_dir
    global marker_vis_dir
    global unique_markers_dir
    global clustree_dir
    global annotation_dir
    global scMCA_annotation_dir
    global fractions_dir
    global gsea_dir
    global scdiffcom_dir
    global cellphonedb_dir
    global gprofiler_dir
    global revigo_dir
    global sc_clustering_dir
    # checkpoint
    checkpoint_dir = output_dir + '/0_checkpoint'
    create_dir(checkpoint_dir)
    # Cellphonedb
    create_dir(cpdb_dir)
    # preprocessing
    preprocessing_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/1_preprocessing'
    create_dir(preprocessing_dir)
    # normalization
    normalization_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/2_normalization'
    create_dir(normalization_dir)
    # feature selection
    feature_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/3_feature_selection'
    create_dir(feature_dir)
    # dimensionality reduction
    dim_reduction_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/4_dim_reduction'
    create_dir(dim_reduction_dir)
    # clustering
    clustering_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/5_embedding_clustering'
    create_dir(clustering_dir)
    # confouding factors
    confounding_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/6_confounding_factors'
    create_dir(confounding_dir)
    # lineage doublets
    lineage_doublets_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_lineage_doublets'
    create_dir(lineage_doublets_dir)
    # lineage marker genes
    lineage_marker_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_lineage_markers'
    create_dir(lineage_marker_dir)
    # putative lineage
    putative_lineage_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_putative_lineage'
    create_dir(putative_lineage_dir)
    # dotplot dir
    dotplot_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_dotplots'
    create_dir(dotplot_dir)
    # scatter plot dir
    scatterplot_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_scatterplots'
    create_dir(scatterplot_dir)
    # violin plot dir
    violinplot_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_violinplots'
    create_dir(violinplot_dir)
    # HCA plot dir
    hca_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_hca'
    create_dir(hca_dir)
    # bar plot dir
    barplot_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_barplots'
    create_dir(barplot_dir)
    # venn diagram dir
    venn_diagram_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_venn_diagrams'
    create_dir(venn_diagram_dir)
    # marker genes
    marker_genes_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/7_marker_genes'
    create_dir(marker_genes_dir)
    # marker genes vizualization
    marker_vis_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/8_marker_visualization'
    create_dir(marker_vis_dir)
    # unique markers info
    unique_markers_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/9_unique_markers'
    create_dir(unique_markers_dir)
    # cluster trees
    clustree_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/10_cluster_trees'
    create_dir(clustree_dir)
    # annotations
    annotation_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/11_annotations'
    create_dir(annotation_dir)
    # scMCA annotations
    scMCA_annotation_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/12_scMCA_annotations'
    create_dir(scMCA_annotation_dir)
    # Cells and sample fractions
    fractions_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/13_cell_sample_fractions'
    create_dir(fractions_dir)
    # GSEA tests
    gsea_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/14_gsea'
    create_dir(gsea_dir)
    # scDiffCom results
    scdiffcom_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/15_scDiffCom'
    create_dir(scdiffcom_dir)
    # cellphonedb results
    cellphonedb_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/16_CellPhoneDB'
    create_dir(cellphonedb_dir)
    # gprofiler results
    gprofiler_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/17_gProfiler'
    create_dir(gprofiler_dir)
    # revigo results
    revigo_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/18_revigo'
    create_dir(revigo_dir)
    create_dir(f"{revigo_dir}/input")
    # SC clustering data
    sc_clustering_dir = output_dir + f'/{cell_suffix}_{f_suffix}' + '/19_SC_embedding_clustering'
    create_dir(sc_clustering_dir)
    

def create_dir(dest: str) -> None:
    import os
    if not os.path.exists(dest):
        os.makedirs(dest)
        print("Created dir: " + dest)


init()
