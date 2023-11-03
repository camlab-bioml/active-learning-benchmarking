def get_maj_min_cell_types(modality, similarity):
    if modality == 'CyTOF':
        if similarity == 'similar': 
            maj = "Intermediate Monocytes"
            minor = "Eosinophils"
        elif similarity == "different": 
            maj = "Intermediate Monocytes"
            minor = "IgDpos IgMpos B cells"
    
    elif modality == "scRNASeq":
        if similarity == 'similar':
            maj = "HCC1937"
            minor = "CAL851"
        elif similarity == "different":
            maj = "HCC1937"
            minor = "MDAMB468"

    elif modality == "snRNASeq":
        if similarity == "similar":
            maj = "Tumor"
            minor = "Atypical_Ductal"
        elif similarity == "different":
            maj = "Tumor"
            minor = "Immune"

    elif modality == "scRNALung":
        if similarity == "similar":
            maj = "HCC827"
            minor = "H1975"
        elif similarity == "different":
            maj = "HCC827"
            minor = "A549"

    elif modality == "liverAtlas":
        if similarity == "similar":
            maj = "T cells"
            minor = "Resident NK"
        elif similarity == "different":
            maj = "T cells"
            minor = "Mono+mono derived cells"
    
    elif modality == "tabulaVasc":
        if similarity == "similar":
            maj = "smooth muscle cell"
            minor = "pericyte cell"
        elif similarity == "different":
            maj = "smooth muscle cell"
            minor = "endothelial cell"
    
    return maj, minor

similarity_vals = ['similar', 'different']
dataset_imbalance = ['balanced', 'imbalanced']

small_selection_strat_dict = {
    'NoMarkerSeurat-clustering': {'initial': ['NA'],
                        'neighbors': Seurat_neighbors,
                        'res': Seurat_resolution,
                        'set': ['NA'],
                        'strategy': 'NA',
                        'AL_alg': 'NA',
                        'random_selection': 'NA',
                        'corruption': [0]},
    'MarkerSeurat-clustering': {'initial': ['NA'],
                        'neighbors': Seurat_neighbors,
                        'res': Seurat_resolution,
                        'set': ['NA'],
                        'strategy': 'NA',
                        'AL_alg': 'NA',
                        'random_selection': 'NA',
                        'corruption': [0]},
    'random': {'initial': ['NA'],
               'neighbors': ['NA'],
               'res': ['NA'],
               'set': random_sets,
               'strategy': 'NA',
               'AL_alg': 'NA',
               'random_selection': 'NA',
               'corruption': [0]},
    'Active-Learning_entropy': {'neighbors': ['NA'],
               'initial': initial_selections,
               'res': ['NA'],
               'set': ['NA'],
               'strategy': ['0.75_quant_entropy', '0.95_quant_entropy', 'highest_entropy'],
               'AL_alg': AL_methods,
               'random_selection': [0],
               'corruption': [0]},
    'Active-Learning_maxp': {'neighbors': ['NA'],
               'initial': initial_selections,
               'res': ['NA'],
               'set': ['NA'],
               'strategy': ['0.05_quant_maxp', '0.25_quant_maxp', 'lowest_maxp'],
               'AL_alg': AL_methods,
               'random_selection': [0],
               'corruption': [0]}
}

scRNA_imb_predictions = []
scRNA_imb = [expand(output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{{modality}}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                        bal = dataset_imbalance,
                        similarity = similarity_vals,
                        initial = small_selection_strat_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = small_selection_strat_dict[select]['strategy'],
                        AL_alg = small_selection_strat_dict[select]['AL_alg'],
                        rand = small_selection_strat_dict[select]['random_selection'],
                        corrupt = small_selection_strat_dict[select]['corruption'],
                        neighbors = small_selection_strat_dict[select]['neighbors'], 
                        res = small_selection_strat_dict[select]['res'], 
                        method = evaluation_methods_dict['scRNASeq'],
                        cell_num = [100]) 
                        for select in small_selection_strat_dict.keys()]
for element in scRNA_imb:
    scRNA_imb_predictions.extend(element)

CyTOF_imb_predictions = []
CyTOF_imb = [expand(output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{{modality}}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                        bal = dataset_imbalance,
                        similarity = similarity_vals,
                        initial = small_selection_strat_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = small_selection_strat_dict[select]['strategy'],
                        AL_alg = small_selection_strat_dict[select]['AL_alg'],
                        rand = small_selection_strat_dict[select]['random_selection'],
                        corrupt = small_selection_strat_dict[select]['corruption'],
                        neighbors = small_selection_strat_dict[select]['neighbors'], 
                        res = small_selection_strat_dict[select]['res'], 
                        method = evaluation_methods_dict['CyTOF'],
                        cell_num = [100]) 
                        for select in small_selection_strat_dict.keys()]

for element in CyTOF_imb:
    CyTOF_imb_predictions.extend(element)

rna_imb_pred_alluvs = []
rna_imb_all = [
    expand(output + 'figures/imbalanced/prediction_alluvials/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf',
        modality = ['scRNASeq', 'snRNASeq', 'CyTOF',],
        bal = dataset_imbalance,
        similarity = similarity_vals,
        initial = small_selection_strat_dict[select]['initial'],
        selection_procedure = [select],
        strat = small_selection_strat_dict[select]['strategy'],
        AL_alg = small_selection_strat_dict[select]['AL_alg'],
        rand = small_selection_strat_dict[select]['random_selection'],
        corruption = small_selection_strat_dict[select]['corruption'],
        knn = small_selection_strat_dict[select]['neighbors'],
        res = small_selection_strat_dict[select]['res']) for select in small_selection_strat_dict.keys()
]
for element in rna_imb_all:
    rna_imb_pred_alluvs.extend(element)

imbalance = {
    'CyTOF': expand(CyTOF_imb_predictions, modality = ['CyTOF'], s = train_test_seeds),
    #'RNA': expand(scRNA_imb_predictions, modality = ['scRNASeq', 'snRNASeq'], s = train_test_seeds),
    'acc': expand(output + 'imbalance/acc/imbalance-acc-{modality}.tsv', modality = 'liverAtlas'),
    #'gt_predicted_alluv': rna_imb_pred_alluvs,
    #'test_freq': output + 'imbalance/cell_type_freqs/test-dataset-cellType-freq.tsv',
    #'report': output + 'reports/imbalance/benchmark-imbalance.html'
}

def get_imbalanced_sce(mod):
    if mod == "snRNASeq":
        sce = 'data/{modality}/{modality}-full-imbalance.rds'
    else:
        sce = 'data/{modality}/{modality}-full.rds'
    return sce

### [ CREATE DATASETS ] #####
rule create_imbalanced_datasets:
    input:
        sce = lambda wildcards: get_imbalanced_sce(wildcards.modality)
    params:
        majority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[0],
        minority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[1]
    output:
        sce = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        rem_sce = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-test-seed-{s}.rds',
        test_set_tsv = 'data/{modality}/imbalance/test-breakdown/{bal}-{similarity}-imbalance-test-seed-{s}.tsv'
    script:
        'imbalance/create-imbalanced-datasets.R'

# rule combine_test_celltype_frequencies:
#     input:
#         tsvs = expand('data/{modality}/imbalance/test-breakdown/{bal}-{similarity}-imbalance-test-seed-{s}.tsv',
#             modality = ['scRNASeq', 'snRNASeq', 'CyTOF'], bal = dataset_imbalance, similarity = similarity_vals, s = train_test_seeds)
#     priority: 2
#     output:
#         tsv = output + 'imbalance/cell_type_freqs/test-dataset-cellType-freq.tsv'
#     script:
#         'imbalance/combine-cellType-freqs.R'

### [ ACTIVE LEARNING CELL SELECTION ] #####
rule Create_active_learning_ground_truth_imbalance:
    input:
        markers = 'markers/{modality}-{similarity}.yml',
        expression = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds'
    priority: 1
    resources: mem_mb=3000
    params:
        max_cell_num = 100
    output:
        assignments = output + 'data/imbalance-{similarity}-{bal}/{modality}/{AL_type}/Init-{initial}-strat-{strat}-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{modality}-knn_neighbors-NA-resolution-NA-iterations_set-full-seed-{s}.tsv',
        entropy = output + 'data/imbalance-{similarity}-{bal}/{modality}/{AL_type}/Init-{initial}-strat-{strat}-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv'
    script:
        'cell-type-assignment/simulate-active-learner.R'

rule create_AL_training_batches_imbalance:
    # This rule subsets the training data by selecting the first n cells
    input:
        assignment = output + 'data/imbalance-{similarity}-{bal}/{modality}/{AL_type}/Init-{initial}-strat-{strat}-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{modality}-knn_neighbors-NA-resolution-NA-iterations_set-full-seed-{s}.tsv'
    output:
        split = output + 'data/imbalance-{similarity}-{bal}/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{modality}-knn_neighbors-NA-resolution-NA-seed-{s}-{subset_val}_cells.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner-n-cells.R'


### [ SEURAT CLUSTERING CELL SELECTION ] #####
rule Seurat_clustering_CyTOF_imbalance:
    input:
        training_rds = 'data/CyTOF/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        markers = 'markers/CyTOF-{similarity}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/{bal}-diagnostics/CyTOF-{bal}-{similarity}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/{bal}-diagnostics/CyTOF-{bal}-{similarity}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "CyTOF"
    output:
        cluster_umap_pdf = output + 'figures/{bal}/CyTOF-{similarity}-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/{bal}/CyTOF-{similarity}-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance-{similarity}-{bal}/CyTOF/cluster-and-interpret/CyTOF-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/{bal}/CyTOF-{similarity}-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_scRNASeq_imbalance:
    input:
        training_rds = 'data/scRNASeq/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        markers = 'markers/scRNASeq-{similarity}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/{bal}-diagnostics/scRNASeq-{bal}-{similarity}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/{bal}-diagnostics/scRNASeq-{bal}-{similarity}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "scRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/{bal}/scRNASeq-{similarity}-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/{bal}/scRNASeq-{similarity}-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance-{similarity}-{bal}/scRNASeq/cluster-and-interpret/scRNASeq-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/{bal}/scRNASeq-{similarity}-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_snRNASeq_imbalance:
    input:
        training_rds = 'data/snRNASeq/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        markers = 'markers/snRNASeq-{similarity}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/{bal}-diagnostics/snRNASeq-{similarity}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/{bal}-diagnostics/snRNASeq-{similarity}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "snRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/{bal}/snRNASeq-{similarity}-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/{bal}/snRNASeq-{similarity}-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance-{similarity}-{bal}/snRNASeq/cluster-and-interpret/snRNASeq-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/{bal}/snRNASeq-{similarity}-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_scRNALung_imbalance:
    input:
        training_rds = 'data/scRNALung/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        markers = 'markers/scRNALung-{similarity}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/{bal}-diagnostics/scRNALung-{similarity}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/{bal}-diagnostics/scRNALung-{similarity}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "scRNALung"
    output:
        cluster_umap_pdf = output + 'figures/{bal}/scRNALung-{similarity}-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/{bal}/scRNALung-{similarity}-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance-{similarity}-{bal}/scRNALung/cluster-and-interpret/scRNALung-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/{bal}/scRNALung-{similarity}-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_liverAtlas_imbalance:
    input:
        training_rds = 'data/liverAtlas/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        markers = 'markers/liverAtlas-{similarity}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/{bal}-diagnostics/liverAtlas-{similarity}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/{bal}-diagnostics/liverAtlas-{similarity}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "liverAtlas"
    output:
        cluster_umap_pdf = output + 'figures/{bal}/liverAtlas-{similarity}-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/{bal}/liverAtlas-{similarity}-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance-{similarity}-{bal}/liverAtlas/cluster-and-interpret/liverAtlas-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/{bal}/liverAtlas-{similarity}-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_tabulaVasc_imbalance:
    input:
        training_rds = 'data/tabulaVasc/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        markers = 'markers/tabulaVasc-{similarity}.yml'
    params:
        positive_markers_diagnostic = output + 'figures/{bal}-diagnostics/tabulaVasc-{similarity}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/{bal}-diagnostics/tabulaVasc-{similarity}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "tabulaVasc"
    output:
        cluster_umap_pdf = output + 'figures/{bal}/tabulaVasc-{similarity}-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/{bal}/tabulaVasc-{similarity}-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance-{similarity}-{bal}/tabulaVasc/cluster-and-interpret/tabulaVasc-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/{bal}/tabulaVasc-{similarity}-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

def get_marker_no_marker_params_clustering_subset(cM):
    if cM == "NoMarkerSeurat":
        param = "wout_markers"
    elif cM == 'MarkerSeurat':
        param = "w_markers"
    return param

rule create_clustering_subsets_imbalance:
    input:
        seurat = output + 'data/imbalance-{similarity}-{bal}/{modality}/cluster-and-interpret/{modality}-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        sce = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds'
    params:
        selection_type = lambda wildcards: get_marker_no_marker_params_clustering_subset(wildcards.clusteringMarkers)
    output:
        ground_truth = output + 'data/imbalance-{similarity}-{bal}/{modality}/{clusteringMarkers}-clustering/Init_NA-strat-NA-rand_sel-NA-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.tsv',
        sce = output + 'data/imbalance-{similarity}-{bal}/{modality}/{clusteringMarkers}-clustering/Init_NA-strat-NA-rand_sel-NA-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.rds'
    script:
        'process-data/select-cluster-subset.R'


### [ RANDOM CELL SELECTION ] #####
rule create_random_subsets_imbalance:
    input:
        sce = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds'
    output:
        sce_subset1 = output + 'data/imbalance-{similarity}-{bal}/{modality}/random/Init_NA-strat-NA-rand_sel-NA-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.rds',
        gt_subset1 = output + 'data/imbalance-{similarity}-{bal}/{modality}/random/Init_NA-strat-NA-rand_sel-NA-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.tsv'
    script:
        'process-data/select-random-subset.R'


### [ CELL TYPE PREDICTION AND MODEL TRAINING ] #####
rule train_and_predict_scmap_imbalance:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality, True),
        train_data = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        test_data = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-test-seed-{s}.rds'
    resources: mem_mb=3000
    output:
        cluster_predictions = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-cluster-predictions-seed-{s}-{cell_num}-cells.tsv',
        sc_predictions = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-sc-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR_imbalance:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality, True),
        train_data = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        test_data = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-test-seed-{s}.rds'
    resources: mem_mb=3000
    output:
        predictions = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleR-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule process_data_for_random_forest_imbalance_train:
    input:
        expression_sce = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds'
    resources:
        mem_mb=5000
    output:
        expression_df = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-expression-df-train-seed-{s}.tsv'
    script:
        'process-data/Create-expression-df.R'

rule process_data_for_random_forest_imbalance_test:
    input:
        expression_sce = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-test-seed-{s}.rds'
    resources:
        mem_mb=5000
    output:
        expression_df = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-expression-df-test-seed-{s}.tsv'
    script:
        'process-data/Create-expression-df.R'

rule train_random_forest_imbalance:
    input:
        train = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-expression-df-train-seed-{s}.tsv',
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality, True)
    output:
        model = output + 'models/{bal}-{similarity}/Init-{initial}-random-forest-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=20000
    log:
        output + 'logs/{bal}-{similarity}-cell-type-predictions/Init-{initial}-random-forest-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule predict_random_forest_imbalance:
    input:
        test = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/{bal}-{similarity}/Init-{initial}-random-forest-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=5000
    output:
        predictions = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-Random-Forest-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule train_SVM_imbalance:
    input:
        train = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-expression-df-train-seed-{s}.tsv',
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality, True),
        markers = 'markers/{modality}.yml'
    params:
        majority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[0],
        minority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[1]
    output:
        model = output + 'models/{bal}-{similarity}/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=5000
    log:
        output + 'logs/{bal}-{similarity}-cell-type-predictions/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/SVM-train.py'

rule predict_SVM_imbalance:
    input:
        test = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/{bal}-{similarity}/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    params:
        majority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[0],
        minority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[1]
    resources:
        mem_mb=2000
    output:
        predictions = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-SVM-rejection-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-SVM.py'

rule train_and_predict_SCN_imbalance:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality, True),
        train_data = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        test_data = 'data/{modality}/imbalance/{similarity}/{bal}-{similarity}-imbalance-test-seed-{s}.rds'
    output:
        predictions = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleCellNet-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleCellNet.R'


rule CyTOF_LDA_imbalance:
    input:
        training_rds = 'data/CyTOF/imbalance/{similarity}/{bal}-{similarity}-imbalance-train-seed-{s}.rds',
        annotation_rds = 'data/CyTOF/imbalance/{similarity}/{bal}-{similarity}-imbalance-test-seed-{s}.rds',
        labels = lambda wildcards: get_labels(wildcards.selection_procedure, 'CyTOF', True),
    resources: mem_mb=3000
    output:
        prediction = output + 'imbalance/rare-subtype-benchmarking/{bal}-{similarity}/Init-{initial}-CyTOF-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-CyTOF-LDA-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'

rule imbalance_compare_gt_to_predicted:
    input:
        sce = lambda wildcards: get_imbalanced_sce(wildcards.modality),
        predicted = lambda wildcards: get_gt_predicted_input_files(wildcards.modality, False, True)
    output:
        pdf = output + 'figures/imbalanced/prediction_alluvials/{bal}-{similarity}/Init-{initial}-{modality}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-gt-against-predicted.R'

def get_imbalanced_predictions(mod):
    if mod == "CyTOF":
        files = CyTOF_imb_predictions
    else:
        files = scRNA_imb_predictions
    
    return files

rule calculate_acc:
    input:
        predictions = lambda wildcards: get_imbalanced_predictions(wildcards.modality),
        sce = lambda wildcards: get_imbalanced_sce(wildcards.modality)
    params:
        majority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[0],
        minority = lambda wildcards: get_maj_min_cell_types(wildcards.modality, wildcards.similarity)[1],
    output:
        acc = output + 'imbalance/acc/imbalance-acc-{modality}-{similarity}-seed-{s}.tsv'
    script:
        'imbalance/calculate-acc.R'

rule combine_accs:
    input:
        accs = expand(output + 'imbalance/acc/imbalance-acc-{{modality}}-{similarity}-seed-{s}.tsv', similarity = similarity_vals, s = train_test_seeds)
    output:
        acc = output + 'imbalance/acc/imbalance-acc-{modality}.tsv'
    script:
        'imbalance/combine-accs.R'


rule benchmark_imbalance:   
    input:
        sc_acc = output + 'imbalance/acc/imbalance-acc-scRNASeq.tsv',
        sn_acc = output + 'imbalance/acc/imbalance-acc-snRNASeq.tsv',
        cy_acc = output + 'imbalance/acc/imbalance-acc-CyTOF.tsv',
        sc_lung = output + 'imbalance/acc/imbalance-acc-scRNALung.tsv',
    params:
        output_dir = output + 'reports/imbalance/'
    output:
        html = output + 'reports/imbalance/benchmark-imbalance.html'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/imbalance/evaluation.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(sn_acc = '{input.sn_acc}', sc_acc = '{input.sc_acc}', cy_acc = '{input.cy_acc}', lung_acc = '{input.sc_lung}'))\" "
