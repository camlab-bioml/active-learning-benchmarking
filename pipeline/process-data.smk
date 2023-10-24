
celllines = ["AU565", "HCC1937", "HCC38", "MDAMB468", "EFM19", "HCC1187", 'JIMT1', "MDAMB361", 'HCC1500', "HCC70", "CAL51", "MDAMB415", 'BT549', "BT483", "MDAMB436", "DU4475", "HS578T", "MCF7", "CAMA1", "BT20", "T47D", "EVSAT", "HDQP1", "BT474", "CAL851", "HCC1143", "MCF12A", "HCC1954", "KPL1", "ZR751", "MX1", "MDAMB453"]

process_data_output = {
    'train_test_split': expand('data/{modality}/{modality}-{split}-seed-{s}.rds', modality = modalities, split = data_splits, s = train_test_seeds),
    'scRNA_expression_df': expand('data/{modality}/{modality}-expression-df-{split}-seed-{s}.tsv', modality = modalities, split = data_splits, s = train_test_seeds),
    'celllines': 'data/scRNASeq/scRNASeq-cellLines-full.rds',
    # 'markers': expand(output + 'figures/{cell_line_mod}-markers/{plot}.pdf', cell_line_mod = ['scRNASeq', 'scRNALung'], plot = ['heatmap', 'lfcs-dot-plot']),
    
    # Random subsets
    'sce_subset1': expand(output + 'data/{modality}/random/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.rds', modality = modalities, cell_num = cell_numbers, s = train_test_seeds),
    'gt_subset1': expand(output + 'data/{modality}/random/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.tsv', modality = modalities, cell_num = cell_numbers, s = train_test_seeds),
    
    # Seurat clustering subsets
    'seu_sce': expand(output + 'data/{modality}/{clusteringMarkers}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.rds', 
        modality = modalities, clusteringMarkers = ['NoMarkerSeurat-clustering', 'MarkerSeurat-clustering'], neighbors = Seurat_neighbors, res = Seurat_resolution, cell_num = cell_numbers, s = train_test_seeds),
    'gt_seu': expand(output + 'data/{modality}/{clusteringMarkers}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.tsv', 
        modality = modalities, clusteringMarkers = ['NoMarkerSeurat-clustering', 'MarkerSeurat-clustering'], neighbors = Seurat_neighbors, res = Seurat_resolution, cell_num = cell_numbers, s = train_test_seeds),
}

# This set of rules is required to process the cell line data and the outputs were used to create a marker file
# rule process_cell_line_mats:
#     input:
#         mat = "data/scRNASeq/matrix.mtx.gz",
#         barcodes = "data/scRNASeq/barcodes.tsv",
#         features = "data/scRNASeq/features.tsv.gz"
#     resources: mem_mb=50000
#     output:
#         sce = temp('data/scRNASeq/{cell_line}-sce.rds')
#     script:
#         'process-data/process-cell-line-data.R'

# rule combine_cell_lines_into_sce:
#     input:
#         mat = expand('data/scRNASeq/{cell_line}-sce.rds', cell_line = celllines)
#     resources: mem_mb=50000
#     output:
#         detected = output + 'figures/scRNASeq-cell-lines-detected.pdf',
#         dim_red = output + 'figures/scRNASeq-cell-lines-dimred.png',
#         sce = 'data/scRNASeq/scRNASeq-cellLines-full.rds'
#     script:
#         'process-data/create-sce.R'

# rule subset_cellLines:
#     input:
#         sce = "data/scRNASeq/scRNASeq-cellLines-full.rds"
#     output:
#         sce = "data/scRNASeq/scRNASeq-full.rds"
#     script:
#         'process-data/subset-celllines.R'

# rule get_markers:
#     input:
#         sce = "data/{cell_line_mod}/{cell_line_mod}-full.rds"
#     output:
#         heatmap = output + 'figures/{cell_line_mod}-markers/heatmap.pdf',
#         lfc_point = output + 'figures/{cell_line_mod}-markers/lfcs-dot-plot.pdf'
#     script:
#         'process-data/scRNASeq-find-markers.R'

rule split_datasets:
    input:
        rds = 'data/{modality}/{modality}-full.rds'
    resources:
        mem_mb=5000
    output:
        train = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test = 'data/{modality}/{modality}-test-seed-{s}.rds',
    script:
        'process-data/split-train-test.R'

rule process_data_for_random_forest:
    input:
        expression_sce = 'data/{modality}/{modality}-{split}-seed-{s}.rds'
    resources:
        mem_mb=5000
    output:
        expression_df = 'data/{modality}/{modality}-expression-df-{split}-seed-{s}.tsv'
    script:
        'process-data/Create-expression-df.R'

rule create_random_subsets:
    input:
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    output:
        sce_subset1 = output + 'data/{modality}/random/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.rds',
        gt_subset1 = output + 'data/{modality}/random/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.tsv',
    script:
        'process-data/select-random-subset.R'

rule create_seurat_clustering_subsets_no_markers:
    input:
        seurat = output + 'cluster-and-interpret/{modality}/{modality}-NoMarkerSeurat-clustering-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    params:
        selection_type = "wout_markers"
    output:
        ground_truth = output + 'data/{modality}/NoMarkerSeurat-clustering/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.tsv',
        sce = output + 'data/{modality}/NoMarkerSeurat-clustering/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.rds'
    script:
        'process-data/select-cluster-subset.R'

rule create_seurat_clustering_subsets_markers:
    input:
        seurat = output + 'cluster-and-interpret/{modality}/{modality}-MarkerSeurat-clustering-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    params:
        selection_type = "w_markers"
    output:
        ground_truth = output + 'data/{modality}/MarkerSeurat-clustering/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.tsv',
        sce = output + 'data/{modality}/MarkerSeurat-clustering/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.rds'
    script:
        'process-data/select-cluster-subset.R'
