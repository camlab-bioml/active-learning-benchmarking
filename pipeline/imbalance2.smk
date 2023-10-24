
CyTOF_imb2_predictions = []
CyTOF_imb2 = [expand(output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-CyTOF-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                        bal = dataset_imbalance,
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

for element in CyTOF_imb2:
    CyTOF_imb2_predictions.extend(element)

scRNA_imb_predictions2 = []
scRNA_imb2 = [expand(output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{{modality}}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                        bal = dataset_imbalance,
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
for element in scRNA_imb2:
    scRNA_imb_predictions2.extend(element)

imbalance2 = {
    #'imbalance2_datasets': expand('data/{modality}/imbalance2/imbalanced-train-seed-{s}.rds', modality = ['CyTOF', 'snRNASeq', 'tabulaVasc'], s = train_test_seeds),
    'CyTOF': expand(CyTOF_imb2_predictions, s = train_test_seeds),
    'scRNASeq': expand(scRNA_imb_predictions2, modality = 'scRNASeq', s = train_test_seeds),
    'snRNASeq': expand(scRNA_imb_predictions2, modality = 'snRNASeq', s = train_test_seeds),
    'scRNALung': expand(scRNA_imb_predictions2, modality = 'scRNALung', s = train_test_seeds),
    'liverAtlas': expand(scRNA_imb_predictions2, modality = 'liverAtlas', s = train_test_seeds),
    'tabulaVasc': expand(scRNA_imb_predictions2, modality = 'tabulaVasc', s = train_test_seeds),
    'imb2_acc': expand(output + 'imbalance2/acc/imbalance-acc-{modality}.tsv', modality = modalities)
}

### Remove these functions after
def get_labels_imbalance2(procedure, mod):
    if procedure == "Active-Learning_entropy" or procedure == "Active-Learning_maxp":
        path = output + 'data/imbalance2-{{bal}}/{modality}/{{selection_procedure}}/AL-batches-subset/Init-{{initial}}-strat-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-{modality}-knn_neighbors-{{neighbors}}-resolution-{{res}}-seed-{{s}}-{{cell_num}}_cells.tsv'
    else:
        path = output + 'data/imbalance2-{{bal}}/{modality}/{{selection_procedure}}/Init_{{initial}}-strat-{{strat}}-rand_sel-{{rand}}-corr-{{corrupt}}-knn_neighbors-{{neighbors}}-resolution-{{res}}-seed-{{s}}-{{cell_num}}_cells.tsv'
    path = expand(path, modality = mod)

    return path
### End of functions to remove


def get_second_imb_cell_types(modality):
    if modality == 'CyTOF':
        cell_types = ['IgDpos IgMpos B cells', 'Classical Monocytes', 'Intermediate Monocytes', 'Eosinophils', 'IgD- IgMpos B cells']
        maj_cell_type = 'IgDpos IgMpos B cells'
    elif modality == 'scRNASeq':
        cell_types = ['MDAMB468', 'CAL51', 'HCC1937', 'CAL851', 'MCF7']
        maj_cell_type = 'MDAMB468'
    elif modality == 'snRNASeq':
        cell_types = ['Tumor', 'Fibroblast', 'Immune', 'Ductal', 'Endothelial']
        maj_cell_type = 'Tumor'
    elif modality == 'scRNALung':
        cell_types = ['A549', 'H838', 'H2228', 'HCC827', 'H1975']
        maj_cell_type = 'A549'
    elif modality == 'liverAtlas':
        cell_types = ['T cells', 'Resident NK', 'Mono+mono derived cells', 'Macrophages', 'Neutrophils']
        maj_cell_type = 'T cells'
    elif modality == 'tabulaVasc':
        cell_types = ['fibroblast', 'macrophage', 'smooth muscle cell', 'endothelial cell', 'pericyte cell']
        maj_cell_type = 'fibroblast'
    
    return cell_types, maj_cell_type

def imbal2_mem(mod, base, mult = 6):
    if mod == 'tabulaVasc':
        base = base * mult
    elif mod == 'liverAtlas':
        base = base * 5
    return base

rule create_2nd_imbalanced_datasets:
    input:
        sce = 'data/{modality}/{modality}-full.rds'
    params:
        cell_types = lambda wildcards: get_second_imb_cell_types(wildcards.modality)[0],
        maj_cell_type = lambda wildcards: get_second_imb_cell_types(wildcards.modality)[1],
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 500)
    output:
        imbal_sce_train = 'data/{modality}/imbalance2/imbalanced-train-seed-{s}.rds',
        imbal_sce_test = 'data/{modality}/imbalance2/imbalanced-test-seed-{s}.rds',
        bal_sce_train = 'data/{modality}/imbalance2/balanced-train-seed-{s}.rds',
        bal_sce_test = 'data/{modality}/imbalance2/balanced-test-seed-{s}.rds'
    script:
        'imbalance/create-second-imbalanced-datasets.R'


rule Create_active_learning_ground_truth_imbalance2:
    input:
        markers = 'markers/imbalance2/{modality}.yml',
        expression = 'data/{modality}/imbalance2/{bal}-train-seed-{s}.rds'
    priority: 1
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 1600)
    params:
        max_cell_num = 100
    output:
        assignments = output + 'data/imbalance2-{bal}/{modality}/{AL_type}/Init-{initial}-strat-{strat}-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{modality}-knn_neighbors-NA-resolution-NA-iterations_set-full-seed-{s}.tsv',
        entropy = output + 'data/imbalance2-{bal}/{modality}/{AL_type}/Init-{initial}-strat-{strat}-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv'
    script:
        'cell-type-assignment/simulate-active-learner.R'

rule create_AL_training_batches_imbalance2:
    # This rule subsets the training data by selecting the first n cells
    input:
        assignment = output + 'data/imbalance2-{bal}/{modality}/{AL_type}/Init-{initial}-strat-{strat}-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{modality}-knn_neighbors-NA-resolution-NA-iterations_set-full-seed-{s}.tsv',
    output:
        split = output + 'data/imbalance2-{bal}/{modality}/{AL_type}/AL-batches-subset/Init-{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-{modality}-knn_neighbors-NA-resolution-NA-seed-{s}-{subset_val}_cells.tsv'
    script:
        'cell-type-assignment/subset-simulated-active-learner-n-cells.R'


### [ SEURAT CLUSTERING CELL SELECTION ] #####
rule Seurat_clustering_CyTOF_imbalance2:
    input:
        training_rds = 'data/CyTOF/imbalance2/{bal}-train-seed-{s}.rds',
        markers = 'markers/imbalance2/CyTOF.yml'
    resources:
        mem_mb=500
    params:
        positive_markers_diagnostic = output + 'figures/imbalance2-{bal}/CyTOF-{bal}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/imbalance2-{bal}/CyTOF-{bal}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "CyTOF"
    output:
        cluster_umap_pdf = output + 'figures/imbalance2-{bal}/CyTOF-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/imbalance2-{bal}/CyTOF-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance2-{bal}/CyTOF/cluster-and-interpret/CyTOF-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/imbalance2-{bal}/CyTOF-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_scRNASeq_imbalance2:
    input:
        training_rds = 'data/scRNASeq/imbalance2/{bal}-train-seed-{s}.rds',
        markers = 'markers/imbalance2/scRNASeq.yml'
    resources:
        mem_mb=500
    params:
        positive_markers_diagnostic = output + 'figures/imbalance2-{bal}/scRNASeq-{bal}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/imbalance2-{bal}/scRNASeq-{bal}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "scRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/imbalance2-{bal}/scRNASeq-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/imbalance2-{bal}/scRNASeq-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance2-{bal}/scRNASeq/cluster-and-interpret/scRNASeq-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/imbalance2-{bal}/scRNASeq-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_snRNASeq_imbalance2:
    input:
        training_rds = 'data/snRNASeq/imbalance2/{bal}-train-seed-{s}.rds',
        markers = 'markers/imbalance2/snRNASeq.yml'
    resources:
        mem_mb=500
    params:
        positive_markers_diagnostic = output + 'figures/imbalance2-{bal}/snRNASeq-{bal}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/imbalance2-{bal}/snRNASeq-{bal}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "snRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/imbalance2-{bal}/snRNASeq-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/imbalance2-{bal}/snRNASeq-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance2-{bal}/snRNASeq/cluster-and-interpret/snRNASeq-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/imbalance2-{bal}/snRNASeq-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_scRNALung_imbalance2:
    input:
        training_rds = 'data/scRNALung/imbalance2/{bal}-train-seed-{s}.rds',
        markers = 'markers/imbalance2/scRNALung.yml'
    resources:
        mem_mb=500
    params:
        positive_markers_diagnostic = output + 'figures/imbalance2-{bal}/scRNALung-{bal}-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/imbalance2-{bal}/scRNALung-{bal}-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "scRNALung"
    output:
        cluster_umap_pdf = output + 'figures/imbalance2-{bal}/scRNALung-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/imbalance2-{bal}/scRNALung-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance2-{bal}/scRNALung/cluster-and-interpret/scRNALung-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/imbalance2-{bal}/scRNALung-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_liverAtlas_imbalance2:
    input:
        training_rds = 'data/liverAtlas/imbalance2/{bal}-train-seed-{s}.rds',
        markers = 'markers/imbalance2/liverAtlas.yml'
    resources:
        mem_mb=500
    params:
        positive_markers_diagnostic = output + 'figures/imbalance2-{bal}/liverAtlas-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/imbalance2-{bal}/liverAtlas-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "liverAtlas"
    output:
        cluster_umap_pdf = output + 'figures/imbalance2-{bal}/liverAtlas-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/imbalance2-{bal}/liverAtlas-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance2-{bal}/liverAtlas/cluster-and-interpret/liverAtlas-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/imbalance2-{bal}/liverAtlas-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_tabulaVasc_imbalance2:
    input:
        training_rds = 'data/tabulaVasc/imbalance2/{bal}-train-seed-{s}.rds',
        markers = 'markers/imbalance2/tabulaVasc.yml'
    resources:
        mem_mb=5000
    params:
        positive_markers_diagnostic = output + 'figures/imbalance2-{bal}-diagnostics/tabulaVasc-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/imbalance2-{bal}-diagnostics/tabulaVasc-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "tabulaVasc"
    output:
        cluster_umap_pdf = output + 'figures/imbalance2-{bal}/tabulaVasc-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/imbalance2-{bal}/tabulaVasc-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'data/imbalance2-{bal}/tabulaVasc/cluster-and-interpret/tabulaVasc-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/imbalance2-{bal}/tabulaVasc-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

def get_marker_no_marker_params_clustering_subset(cM):
    if cM == "NoMarkerSeurat":
        param = "wout_markers"
    elif cM == 'MarkerSeurat':
        param = "w_markers"
    return param

rule create_clustering_subsets_imbalance2:
    input:
        seurat = output + 'data/imbalance2-{bal}/{modality}/cluster-and-interpret/{modality}-{clusteringMarkers}-assignments-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        sce = 'data/{modality}/imbalance2/{bal}-train-seed-{s}.rds'
    params:
        selection_type = lambda wildcards: get_marker_no_marker_params_clustering_subset(wildcards.clusteringMarkers)
    output:
        ground_truth = output + 'data/imbalance2-{bal}/{modality}/{clusteringMarkers}-clustering/Init_NA-strat-NA-rand_sel-NA-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.tsv',
        sce = output + 'data/imbalance2-{bal}/{modality}/{clusteringMarkers}-clustering/Init_NA-strat-NA-rand_sel-NA-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells.rds'
    script:
        'process-data/select-cluster-subset.R'


### [ RANDOM CELL SELECTION ] #####
rule create_random_subsets_imbalance2:
    input:
        sce = 'data/{modality}/imbalance2/{bal}-train-seed-{s}.rds'
    output:
        sce_subset1 = output + 'data/imbalance2-{bal}/{modality}/random/Init_NA-strat-NA-rand_sel-NA-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.rds',
        gt_subset1 = output + 'data/imbalance2-{bal}/{modality}/random/Init_NA-strat-NA-rand_sel-NA-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-{cell_num}_cells.tsv'
    script:
        'process-data/select-random-subset.R'

### [ CELL TYPE PREDICTION AND MODEL TRAINING ] #####
rule train_and_predict_scmap_imbalance2:
    input:
        annotation = lambda wildcards: get_labels_imbalance2(wildcards.selection_procedure, wildcards.modality),
        train_data = 'data/{modality}/imbalance2/{bal}-train-seed-{s}.rds',
        test_data = 'data/{modality}/imbalance2/{bal}-test-seed-{s}.rds'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 2000, 9)
    output:
        cluster_predictions = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-cluster-predictions-seed-{s}-{cell_num}-cells.tsv',
        sc_predictions = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-sc-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR_imbalance2:
    input:
        annotation = lambda wildcards: get_labels_imbalance2(wildcards.selection_procedure, wildcards.modality),
        train_data = 'data/{modality}/imbalance2/{bal}-train-seed-{s}.rds',
        test_data = 'data/{modality}/imbalance2/{bal}-test-seed-{s}.rds'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 2000)
    output:
        predictions = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleR-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule process_data_for_random_forest_imbalance2:
    input:
        expression_sce = 'data/{modality}/imbalance2/{bal}-{split}-seed-{s}.rds'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 2500, 10)
    output:
        expression_df = 'data/{modality}/imbalance2/{bal}-expression-df-{split}-seed-{s}.tsv'
    script:
        'process-data/Create-expression-df.R'

rule train_random_forest_imbalance2:
    input:
        train = 'data/{modality}/imbalance2/{bal}-expression-df-train-seed-{s}.tsv',
        annotation = lambda wildcards: get_labels_imbalance2(wildcards.selection_procedure, wildcards.modality)
    output:
        model = output + 'models/imbalance2/{bal}/Init-{initial}-random-forest-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 2000)
    log:
        output + 'logs/imbalance2-{bal}-cell-type-predictions/Init-{initial}-random-forest-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule predict_random_forest_imbalance2:
    input:
        test = 'data/{modality}/imbalance2/{bal}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/imbalance2/{bal}/Init-{initial}-random-forest-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 10000, 25)
    output:
        predictions = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-Random-Forest-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule train_SVM_imbalance2:
    input:
        train = 'data/{modality}/imbalance2/{bal}-expression-df-train-seed-{s}.tsv',
        annotation = lambda wildcards: get_labels_imbalance2(wildcards.selection_procedure, wildcards.modality),
        markers = 'markers/imbalance2/{modality}.yml'
    output:
        model = output + 'models/imbalance2/{bal}/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 2000)
    log:
        output + 'logs/imbalance2-{bal}-cell-type-predictions/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/SVM-train.py'

rule predict_SVM_imbalance2:
    input:
        test = 'data/{modality}/imbalance2/{bal}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/imbalance2/{bal}/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 4000, 10)
    output:
        predictions = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-SVM-rejection-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-SVM.py'

rule train_and_predict_SCN_imbalance2:
    input:
        annotation = lambda wildcards: get_labels_imbalance2(wildcards.selection_procedure, wildcards.modality),
        train_data = 'data/{modality}/imbalance2/{bal}-train-seed-{s}.rds',
        test_data = 'data/{modality}/imbalance2/{bal}-test-seed-{s}.rds'
    resources:
        mem_mb=lambda wildcards: imbal2_mem(wildcards.modality, 2000)
    output:
        predictions = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-{modality}-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleCellNet-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleCellNet.R'


rule CyTOF_LDA_imbalance2:
    input:
        training_rds = 'data/CyTOF/imbalance2/{bal}-train-seed-{s}.rds',
        annotation_rds = 'data/CyTOF/imbalance2/{bal}-test-seed-{s}.rds',
        labels = lambda wildcards: get_labels_imbalance2(wildcards.selection_procedure, 'CyTOF'),
    resources: mem_mb=3000
    output:
        prediction = output + 'imbalance2/rare-subtype-benchmarking/{bal}/Init-{initial}-CyTOF-sel-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-CyTOF-LDA-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'

def get_imbalance2_predictions(modality):
    if modality == "CyTOF":
        f = expand(CyTOF_imb2_predictions, modality = 'CyTOF', s = train_test_seeds)
    else:
        f = expand(scRNA_imb_predictions2, modality = modality, s = train_test_seeds)
    return f


rule calculate_acc_imbalance2:
    input:
        predictions = lambda wildcards: get_imbalance2_predictions(wildcards.modality),
        sce = 'data/{modality}/{modality}-full.rds'
    params:
        cell_types = lambda wildcards: get_second_imb_cell_types(wildcards.modality)[0],
    output:
        acc = output + 'imbalance2/acc/imbalance-acc-{modality}.tsv'
    script:
        'imbalance2/calculate-acc.R'
