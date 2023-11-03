
marker_rem_percentages = [0, 0.1, 0.25, 0.5, 0.75, 1]
rna_modalities = ['scRNASeq', 'snRNASeq', 'scRNALung', 'liverAtlas', 'tabulaVasc']

AR_remove_markers = {
    'marker_files': expand(output + 'markers/{modality}-percent-rem-{rem_percentage}.yml', modality = modalities, rem_percentage = marker_rem_percentages),
    'pred': expand(output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv',
        initial = 'NA', modality = rna_modalities, selection_procedure = 'MarkerSeurat-clustering', strat = 'NA', AL_alg = 'NA', rand = 0, corrupt = 0, neighbors = 20, res = 0.8, method = evaluation_methods_dict['scRNASeq'], cell_num = 100, rem_percentage = marker_rem_percentages, s = train_test_seeds),
    'pred_cytof': expand(output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-CyTOF-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv',
        initial = 'NA', selection_procedure = 'MarkerSeurat-clustering', strat = 'NA', AL_alg = 'NA', rand = 0, corrupt = 0, neighbors = 20, res = 0.8, method = evaluation_methods_dict['CyTOF'], cell_num = 100, rem_percentage = marker_rem_percentages, s = train_test_seeds),
    'acc': expand(output + 'results/marker_corruption/overall-{modality}-benchmarking-accuracies.tsv', modality = modalities)
}


rule create_marker_files:
    input:
        sce = 'data/{modality}/{modality}-train-seed-0.rds',
        markers = 'markers/{modality}.yml'
    output:
        markers = output + 'markers/{modality}-percent-rem-{rem_percentage}.yml'
    script:
        'marker-removal/remove-marker-AR.R'

rule AR_rem_markers:
    input:
        training_rds = 'data/{modality}/{modality}-train-seed-{s}.rds',
        markers = output + 'markers/{modality}-percent-rem-{rem_percentage}.yml'
    resources:
        mem_mb=5000
    params:
        positive_markers_diagnostic = output + 'figures/marker_corruption/diagnostics/{modality}-percent-rem-{rem_percentage}-[cell_types]-positive-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        negative_markers_diagnostic = output + 'figures/marker_corruption/diagnostics/{modality}-percent-rem-{rem_percentage}-[cell_types]-negative-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        mod = '{modality}'
    output:
        cluster_umap_pdf = output + 'figures/marker_corruption/{modality}-percent-rem-{rem_percentage}-cluster-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        cell_type_umap_pdf = output + 'figures/marker_corruption/{modality}-percent-rem-{rem_percentage}-cell-assignment-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf',
        assignments = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-full-assignments.tsv',
        ground_truth_umap_pdf = output + 'figures/marker_corruption/{modality}-percent-rem-{rem_percentage}-ground-truth-umap-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}.pdf'
    script:
        'cell-type-assignment/Seurat.R'

rule AR_rem_markers_create_subsets_markers:
    input:
        seurat = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-knn_neighbors-20-resolution-0.8-seed-{s}-full-assignments.tsv',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds',
    params:
        selection_type = "w_markers"
    output:
        ground_truth = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-{cell_num}-cells.tsv',
        sce = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-{cell_num}.rds'
    script:
        'process-data/select-cluster-subset.R'


def get_mem_mb(wildcards, attempt):
    return attempt * 5000 + 3000

rule train_and_predict_scmap_marker_rem:
    input:
        annotation = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-100-cells.tsv',
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=3000
    output:
        cluster_predictions = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-cluster-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv',
        sc_predictions = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-scmap-sc-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule train_and_predict_singleR_marker_rem:
    input:
        annotation = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-100-cells.tsv',
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=3000
    output:
        predictions = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleR-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule train_random_forest_marker_rem:
    input:
        train = 'data/{modality}/{modality}-expression-df-train-seed-{s}.tsv',
        annotation = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-100-cells.tsv'
    output:
        model = output + 'models/marker_corruption/random-forest-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.pkl'
    resources:
        mem_mb=get_mem_mb
    log:
        output + 'logs/marker_corruption/cell-type-predictions/random-forest-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule predict_random_forest_marker_rem:
    input:
        test = 'data/{modality}/{modality}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/marker_corruption/random-forest-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.pkl'
    resources:
        mem_mb=get_mem_mb
    output:
        predictions = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-Random-Forest-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule train_SVM_marker_rem:
    input:
        train = 'data/{modality}/{modality}-expression-df-train-seed-{s}.tsv',
        annotation = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-100-cells.tsv',
        markers = 'markers/{modality}.yml'
    output:
        model = output + 'models/marker_corruption/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.pkl'
    resources:
        mem_mb=get_mem_mb
    log:
        output + 'logs/marker_corruption/cell-type-predictions/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.log'
    script:
        'cell-type-assignment/SVM-train.py'

rule predict_SVM_marker_rem:
    input:
        test = 'data/{modality}/{modality}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/marker_corruption/SVM-rejection-Init_{initial}-{modality}-trained-on-{selection_procedure}-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.pkl',
        markers = 'markers/{modality}.yml'
    resources:
        mem_mb=get_mem_mb
    output:
        predictions = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-SVM-rejection-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv'
    script:
        'cell-type-assignment/predict-SVM.py'

rule train_and_predict_SCN_marker_rem:
    input:
        annotation = output + 'marker_corruption/AR/{modality}-percent-rem-{rem_percentage}-seed-{s}-100-cells.tsv',
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    output:
        predictions = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-singleCellNet-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv'
    script:
        'cell-type-assignment/singleCellNet.R'

rule CyTOF_LDA_marker_rem:
    input:
        training_rds = 'data/CyTOF/CyTOF-train-seed-{s}.rds',
        annotation_rds = 'data/CyTOF/CyTOF-test-seed-{s}.rds',
        labels = output + 'marker_corruption/AR/CyTOF-percent-rem-{rem_percentage}-seed-{s}-100-cells.tsv'
    output:
        prediction = output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-CyTOF-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-CyTOF-LDA-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'

def get_AR_rem_predictions(mod):
    if mod == "CyTOF":
        files = expand(output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv',
            modality = mod, initial = 'NA', selection_procedure = 'MarkerSeurat-clustering', strat = 'NA', AL_alg = 'NA', rand = 0, corrupt = 0, neighbors = 20, res = 0.8, method = evaluation_methods_dict['CyTOF'], cell_num = 100, rem_percentage = marker_rem_percentages, s = train_test_seeds)
    else:
        files = expand(output + 'marker_corruption/rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{s}-{cell_num}-cells-percent-rem-{rem_percentage}.tsv',
            modality = mod, initial = 'NA', selection_procedure = 'MarkerSeurat-clustering', strat = 'NA', AL_alg = 'NA', rand = 0, corrupt = 0, neighbors = 20, res = 0.8, method = evaluation_methods_dict['scRNASeq'], cell_num = 100, rem_percentage = marker_rem_percentages, s = train_test_seeds)

    return files


rule overall_benchmark_marker_rem:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        predictions = lambda wildcards: get_AR_rem_predictions(wildcards.modality)
    resources:
        mem_mb=5000
    log:
        output + 'logs/benchmark-predictive-labeling-{modality}.log'
    output:
        acc = output + 'results/marker_corruption/overall-{modality}-benchmarking-accuracies.tsv'
    script:
        'benchmarking/benchmark-AR-marker-rem.R'
