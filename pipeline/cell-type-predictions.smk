
def expand_predictions_by_mod(mod,
                             selection_expansion_dict = selection_expansion_dict,
                             evaluation_methods_dict = evaluation_methods_dict,
                             train_test_seeds = train_test_seeds,
                             cell_numbers = cell_numbers):
    rand_0 = []
    pred_rand_0 = [expand(output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                        modality = mod,
                        initial = selection_expansion_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = selection_expansion_dict[select]['strategy'],
                        AL_alg = selection_expansion_dict[select]['AL_alg'],
                        rand = 0,
                        corrupt = selection_expansion_dict[select]['corruption'],
                        neighbors = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        method = evaluation_methods_dict[mod],
                        #s = train_test_seeds,
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]
    
    for i in pred_rand_0:
        rand_0.extend(i)

    corr_0 = []
    pred_corr_0 = [expand(output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                        modality = mod,
                        initial = selection_expansion_dict[select]['initial'],
                        selection_procedure = [select], 
                        strat = selection_expansion_dict[select]['strategy'],
                        AL_alg = selection_expansion_dict[select]['AL_alg'],
                        rand = selection_expansion_dict[select]['random_selection'],
                        corrupt = 0,
                        neighbors = selection_expansion_dict[select]['neighbors'], 
                        res = selection_expansion_dict[select]['res'], 
                        method = evaluation_methods_dict[mod],
                        #s = train_test_seeds,
                        cell_num = cell_numbers) 
                        for select in selection_expansion_dict.keys()]
    
    for i in pred_corr_0:
        corr_0.extend(i)

    return rand_0 + corr_0

def get_labels(procedure, mod, imbalanced = False, doublet = False):
    if imbalanced:
        if procedure == "Active-Learning_entropy" or procedure == "Active-Learning_maxp":
            path = output + 'data/imbalance-{{similarity}}-{{bal}}/{modality}/{{selection_procedure}}/AL-batches-subset/Init-{{initial}}-strat-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-{modality}-knn_neighbors-{{neighbors}}-resolution-{{res}}-seed-{{s}}-{{cell_num}}_cells.tsv'
        else:
            path = output + 'data/imbalance-{{similarity}}-{{bal}}/{modality}/{{selection_procedure}}/Init_{{initial}}-strat-{{strat}}-rand_sel-{{rand}}-corr-{{corrupt}}-knn_neighbors-{{neighbors}}-resolution-{{res}}-seed-{{s}}-{{cell_num}}_cells.tsv'
        path = expand(path, modality = mod)
    else:
        if procedure == "Active-Learning_entropy" or procedure == "Active-Learning_maxp":
            path = output + 'data/{modality}/{{selection_procedure}}/AL-batches-subset/Init_{{initial}}-strat-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-knn_neighbors-{{neighbors}}-resolution-{{res}}-seed-{{s}}-{{cell_num}}_cells.tsv'
        else:
            path = output + 'data/{modality}/{{selection_procedure}}/Init_{{initial}}-strat-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-knn_neighbors-{{neighbors}}-resolution-{{res}}-seed-{{s}}-{{cell_num}}_cells.tsv'
        
        path = expand(path, modality = mod)
    return path

cell_type_predictions = {
    'scRNASeq': expand_predictions_by_mod("scRNASeq"),
    'snRNASeq': expand_predictions_by_mod("snRNASeq"),
    'CyTOF': expand_predictions_by_mod("CyTOF"),
    'scRNALung': expand_predictions_by_mod("scRNALung"),
    'liverAtlas': expand_predictions_by_mod("liverAtlas"),
    'tabulaVasc': expand_predictions_by_mod("tabulaVasc")
}

def get_mem_mb(wildcards, attempt):
    return attempt * 5000 + 3000

rule train_and_predict_scmap:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality),
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=3000
    output:
        cluster_predictions = output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-cluster-predictions-seed-{s}-{cell_num}-cells.tsv',
        sc_predictions = output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-sc-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality),
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=3000
    output:
        predictions = output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleR-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule train_random_forest:
    input:
        train = 'data/{modality}/{modality}-expression-df-train-seed-{s}.tsv',
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality)
    output:
        model = output + 'models/random-forest-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=get_mem_mb
    log:
        output + 'logs/cell-type-predictions/random-forest-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule predict_random_forest:
    input:
        test = 'data/{modality}/{modality}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/random-forest-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=get_mem_mb
    output:
        predictions = output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-Random-Forest-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule train_SVM:
    input:
        train = 'data/{modality}/{modality}-expression-df-train-seed-{s}.tsv',
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality),
        markers = 'markers/{modality}.yml'
    output:
        model = output + 'models/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl'
    resources:
        mem_mb=get_mem_mb
    log:
        output + 'logs/cell-type-predictions/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.log'
    script:
        'cell-type-assignment/SVM-train.py'

rule predict_SVM:
    input:
        test = 'data/{modality}/{modality}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells.pkl',
        markers = 'markers/{modality}.yml'
    resources:
        mem_mb=get_mem_mb
    output:
        predictions = output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-SVM-rejection-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/predict-SVM.py'

rule train_and_predict_SCN:
    input:
        annotation = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality),
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    output:
        predictions = output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleCellNet-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/singleCellNet.R'

rule Seurat_clustering_scRNASeq:
    input:
        training_rds = 'data/scRNASeq/scRNASeq-train-seed-{s}.rds',
        markers = 'markers/scRNASeq.yml'
    resources:
        mem_mb=5000
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/scRNASeq-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/scRNASeq-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "scRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/scRNASeq-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/scRNASeq-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'cluster-and-interpret/scRNASeq/scRNASeq-{clusteringMarkers}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/scRNASeq-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_snRNASeq:
    input:
        training_rds = 'data/snRNASeq/snRNASeq-train-seed-{s}.rds',
        markers = 'markers/snRNASeq.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/snRNASeq-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/snRNASeq-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "snRNASeq"
    output:
        cluster_umap_pdf = output + 'figures/snRNASeq-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/snRNASeq-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'cluster-and-interpret/snRNASeq/snRNASeq-{clusteringMarkers}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/snRNASeq-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_scRNALung:
    input:
        training_rds = 'data/scRNALung/scRNALung-train-seed-{s}.rds',
        markers = 'markers/scRNALung.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/scRNALung-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/scRNALung-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "scRNALung"
    output:
        cluster_umap_pdf = output + 'figures/scRNALung-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/scRNALung-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'cluster-and-interpret/scRNALung/scRNALung-{clusteringMarkers}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/scRNALung-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_liverAtlas:
    input:
        training_rds = 'data/liverAtlas/liverAtlas-train-seed-{s}.rds',
        markers = 'markers/liverAtlas.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/liverAtlas-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/liverAtlas-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "liverAtlas"
    resources:
        mem_mb=4000
    output:
        cluster_umap_pdf = output + 'figures/liverAtlas-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/liverAtlas-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'cluster-and-interpret/liverAtlas/liverAtlas-{clusteringMarkers}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/liverAtlas-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_tabulaVasc:
    input:
        training_rds = 'data/tabulaVasc/tabulaVasc-train-seed-{s}.rds',
        markers = 'markers/tabulaVasc.yml'
    resources:
        mem_mb=10000
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/tabulaVasc-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/tabulaVasc-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "tabulaVasc"
    output:
        cluster_umap_pdf = output + 'figures/tabulaVasc-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/tabulaVasc-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'cluster-and-interpret/tabulaVasc/tabulaVasc-{clusteringMarkers}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/tabulaVasc-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule Seurat_clustering_CyTOF:
    input:
        training_rds = 'data/CyTOF/CyTOF-train-seed-{s}.rds',
        markers = 'markers/CyTOF.yml'
    params:
        positive_markers_diagnostic = output + 'figures/diagnostics/CyTOF-{clusteringMarkers}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/diagnostics/CyTOF-{clusteringMarkers}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = "CyTOF"
    output:
        cluster_umap_pdf = output + 'figures/CyTOF-{clusteringMarkers}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/CyTOF-{clusteringMarkers}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'cluster-and-interpret/CyTOF/CyTOF-{clusteringMarkers}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.tsv',
        ground_truth_umap_pdf = output + 'figures/CyTOF-{clusteringMarkers}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'


rule CyTOF_LDA:
    input:
        training_rds = 'data/CyTOF/CyTOF-train-seed-{s}.rds',
        annotation_rds = 'data/CyTOF/CyTOF-test-seed-{s}.rds',
        labels = lambda wildcards: get_labels(wildcards.selection_procedure, 'CyTOF'),
    output:
        prediction = output + 'rare-subtype-benchmarking/Init_{initial}-CyTOF-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-CyTOF-LDA-predictions-seed-{s}-{cell_num}-cells.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'
        
