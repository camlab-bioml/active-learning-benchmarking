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


def get_mem_mb(wildcards, attempt):
    return attempt * 5000 + 10000

def get_mem_mb_pred_labs(wildcards, attempt):
    return attempt * 7000 + 5000

def get_mem_mb_classifiers(wildcards, attempt):
    return attempt * 4000 + 1000

def get_pred_labeling(modality, cell_selection = ['top10', 'top50', 'top100'], selection_expansion_dict = selection_expansion_dict, evaluation_methods_dict = evaluation_methods_dict, cell_numbers = cell_numbers):    
    if modality == 'scRNASeq' or modality == 'snRNASeq' or modality == 'scRNALung' or modality == 'liverAtlas' or modality == 'tabulaVasc':
        data = []
        scRNA = [expand(output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{{s}}-method-{method}-predictions-{cell_num}-cells-{cell_selection}.tsv',
                                modality = modality,
                                initial = selection_expansion_dict[select]['initial'],
                                selection_procedure = [select], 
                                strat = selection_expansion_dict[select]['strategy'],
                                AL_alg = selection_expansion_dict[select]['AL_alg'],
                                pred_alg = ['multinom', 'rf'],
                                rand = [0],
                                corrupt = [0],
                                neighbors = selection_expansion_dict[select]['neighbors'], 
                                res = selection_expansion_dict[select]['res'], 
                                method = evaluation_methods_dict[modality],
                                cell_num = cell_numbers,
                                cell_selection = cell_selection) 
                                for select in selection_expansion_dict.keys()]
        for element in scRNA:
            data.extend(element)
    elif modality == 'CyTOF':
        data = []
        CyTOF = [expand(output + 'predictive-labeling-benchmarking2/CyTOF/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{{s}}-method-{method}-predictions-{cell_num}-cells-{cell_selection}.tsv',
                                initial = selection_expansion_dict[select]['initial'],
                                selection_procedure = [select], 
                                strat = selection_expansion_dict[select]['strategy'],
                                AL_alg = selection_expansion_dict[select]['AL_alg'],
                                pred_alg = ['multinom', 'rf'],
                                rand = [0],
                                corrupt = [0],
                                neighbors = selection_expansion_dict[select]['neighbors'], 
                                res = selection_expansion_dict[select]['res'], 
                                method = evaluation_methods_dict['CyTOF'],
                                cell_num = cell_numbers,
                                cell_selection = cell_selection) 
                                for select in selection_expansion_dict.keys()]
        for element in CyTOF:
            data.extend(element)
    return data

# This function gets all the predictions for the baseline
def pred_lab_expand_normal_predictions_by_mod(mod,
                             method_sel,
                             selection_expansion_dict = selection_expansion_dict,
                             evaluation_methods_dict = evaluation_methods_dict,
                             cell_numbers = cell_numbers):
    f = expand(output + 'rare-subtype-benchmarking/Init_{initial}-{modality}-sel_{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-0-corr-0-knn_neighbors-{neighbors}-resolution-{res}-{method}-predictions-seed-{{s}}-{cell_num}-cells.tsv',
                    modality = mod,
                    initial = selection_expansion_dict[method_sel]['initial'],
                    selection_procedure = [method_sel], 
                    strat = selection_expansion_dict[method_sel]['strategy'],
                    AL_alg = selection_expansion_dict[method_sel]['AL_alg'],
                    neighbors = selection_expansion_dict[method_sel]['neighbors'], 
                    res = selection_expansion_dict[method_sel]['res'], 
                    method = evaluation_methods_dict[mod],
                    cell_num = cell_numbers)

    return f


pred_lab2 = {
    #'pred_labs1': expand(get_pred_labeling("CyTOF") + get_pred_labeling("scRNASeq") + get_pred_labeling("snRNASeq"), s = train_test_seeds)# + 
    #'pred_labs2': expand(get_pred_labeling("scRNALung"), s = train_test_seeds),
    'pred_labsliv': expand(get_pred_labeling("liverAtlas"), s = train_test_seeds),
    #'pred_labs_vasc': expand(get_pred_labeling("tabulaVasc"), s = train_test_seeds),
    # 'acc': expand(output + 'new/benchmark-predictive-labeling-{modality}-seed-{s}.tsv', modality = modalities, s = train_test_seeds),
    # 'benchmark': expand(output + 'reports/pred2/benchmark-predictive-labeling-{modality}.html', modality = modalities),
    'benchmark2': expand(output + 'new/pred2/benchmark-predictive-labeling-{modality}.tsv', modality = modalities),
    'acc_pred_lab': output + 'new/pred-labeling-accuracy.tsv',
    'detect_mislab': expand(output + 'identify_mislabelled/{modality}/Pred_alg-{pred_alg}-corr-10-seed-{seed}-cellNum-250-cells.tsv', modality = modalities, pred_alg = ['rf', 'multinom'], seed = train_test_seeds),
}

rule create_predictive_labels:
    input: 
        markers = 'markers/{modality}.yml',
        sce = "data/{modality}/{modality}-train-seed-{s}.rds",
        selected_cells = lambda wildcards: get_labels(wildcards.selection_procedure, wildcards.modality)
    resources:
        mem_mb=get_mem_mb_pred_labs
    output:
        full_pred = output + 'predictive_labs_data2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-full_pred.tsv'
    script:
        'predictive-labeling/select-and-predict-cells2.R'

rule create_percentile_pred_lab_datasets:
    input:
        full_pred = output + 'predictive_labs_data2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-full_pred.tsv'
    output:
        data = output + 'predictive_labs_data_subset/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv'
    script:
        'predictive-labeling/select-cells.R'

rule calculate_pred_lab_acc:
    input:
        sce = "data/{mod}/{mod}-full.rds",
        AL_ent = expand(output + 'predictive_labs_data_subset/{{mod}}/Init_{initial}-sel-Active-Learning_entropy-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_lab}-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{{s}}-{cell_num}_cells-{sel}.tsv',
            initial = initial_selections, strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], AL_alg = AL_methods, pred_lab = ['multinom', 'rf'], cell_num = cell_numbers, sel = ['top10', 'top50', 'top100']),
        AL_maxp = expand(output + 'predictive_labs_data_subset/{{mod}}/Init_{initial}-sel-Active-Learning_maxp-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_lab}-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{{s}}-{cell_num}_cells-{sel}.tsv',
            initial = initial_selections, strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], AL_alg = AL_methods, pred_lab = ['multinom', 'rf'], cell_num = cell_numbers, sel = ['top10', 'top50', 'top100']),
        random = expand(output + 'predictive_labs_data_subset/{{mod}}/Init_NA-sel-random-strat-NA-ALAlg-NA-pred_alg-{pred_lab}-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{{s}}-{cell_num}_cells-{sel}.tsv',
            cell_num = cell_numbers, sel = ['top10', 'top50', 'top100'], pred_lab = ['multinom', 'rf']),
        AR = expand(output + 'predictive_labs_data_subset/{{mod}}/Init_NA-sel-{sel_method}-strat-NA-ALAlg-NA-pred_alg-{pred_lab}-rand_sel-0-corr-0-knn_neighbors-{knn}-resolution-{res}-seed-{{s}}-{cell_num}_cells-{sel}.tsv',
            sel_method = ['NoMarkerSeurat-clustering', 'MarkerSeurat-clustering'], pred_lab = ['multinom', 'rf'], knn = Seurat_neighbors, res = Seurat_resolution, cell_num = cell_numbers, sel = ['top10', 'top50', 'top100']),
    output:
        acc = temp(output + 'new/pred-labeling-{mod}-accuracy-seed-{s}.tsv')
    script:
        'predictive-labeling/calculate-predictive-labeling-accuracy.R'

rule combine_pred_lab_acc:
    input:
        accs = expand(output + 'new/pred-labeling-{modality}-accuracy-seed-{s}.tsv', modality = modalities, s = train_test_seeds)
    output:
        acc = output + 'new/pred-labeling-accuracy.tsv'
    script:
        'benchmarking/combine-pred-labeling-accs.R'

rule pred_lab_train_and_predict_scmap:
    input:
        annotation = output + 'predictive_labs_data_subset/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv',
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=get_mem_mb_classifiers
    output:
        cluster_predictions = output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-scmap-cluster-predictions-{cell_num}-cells-{cell_selection}.tsv',
        sc_predictions = output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-scmap-sc-predictions-{cell_num}-cells-{cell_selection}.tsv'
    script:
        'cell-type-assignment/scmap.R'

rule pred_lab_train_and_predict_singleR:
    input:
        annotation = output + 'predictive_labs_data_subset/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv',
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=get_mem_mb_classifiers
    output:
        predictions = output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-singleR-predictions-{cell_num}-cells-{cell_selection}.tsv'
    script:
        'cell-type-assignment/singleR.R'

rule pred_lab_train_random_forest:
    input:
        train = 'data/{modality}/{modality}-expression-df-train-seed-{s}.tsv',
        annotation = output + 'predictive_labs_data_subset/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv'
    output:
        model = output + 'models_pred_labs2/random-forest-{modality}-trained-on-Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-{cell_selection}.pkl'
    resources:
        mem_mb=20000
    log:
        output + 'logs/predictive-labeling-benchmarking2/random-forest-{modality}-trained-on-Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-{cell_selection}.log'
    script:
        'cell-type-assignment/random-forest-train.py'

rule pred_lab_predict_random_forest:
    input:
        test = 'data/{modality}/{modality}-expression-df-test-seed-{s}.tsv',
        model = output + 'models_pred_labs2/random-forest-{modality}-trained-on-Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-{cell_selection}.pkl'
    resources:
        mem_mb=get_mem_mb_classifiers
    output:
        predictions = output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-Random-Forest-predictions-{cell_num}-cells-{cell_selection}.tsv'
    script:
        'cell-type-assignment/predict-random-forest.py'

rule pred_lab_train_SVM:
    input:
        train = 'data/{modality}/{modality}-expression-df-train-seed-{s}.tsv',
        annotation = output + 'predictive_labs_data_subset/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv',
        markers = 'markers/{modality}.yml'
    resources:
        mem_mb=get_mem_mb_classifiers
    output:
        model = output + 'models/models_pred_labs2/SVM-rejection-{modality}-trained-on-Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-{cell_selection}.pkl'
    script:
        'cell-type-assignment/SVM-train.py'

rule pred_lab_predict_SVM:
    input:
        test = 'data/{modality}/{modality}-expression-df-test-seed-{s}.tsv',
        model = output + 'models/models_pred_labs2/SVM-rejection-{modality}-trained-on-Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}-cells-{cell_selection}.pkl',
        markers = 'markers/{modality}.yml'
    resources:
        mem_mb=get_mem_mb_classifiers
    output:
        predictions = output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-SVM-rejection-predictions-{cell_num}-cells-{cell_selection}.tsv'
    script:
        'cell-type-assignment/predict-SVM.py'

rule pred_lab_train_and_predict_SCN:
    input:
        annotation = output + 'predictive_labs_data_subset/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv',
        train_data = 'data/{modality}/{modality}-train-seed-{s}.rds',
        test_data = 'data/{modality}/{modality}-test-seed-{s}.rds'
    resources:
        mem_mb=get_mem_mb_classifiers
    output:
        predictions = output + 'predictive-labeling-benchmarking2/{modality}/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-singleCellNet-predictions-{cell_num}-cells-{cell_selection}.tsv'
    script:
        'cell-type-assignment/singleCellNet.R'

rule pred_lab_CyTOF_LDA:
    input:
        training_rds = 'data/CyTOF/CyTOF-train-seed-{s}.rds',
        annotation_rds = 'data/CyTOF/CyTOF-test-seed-{s}.rds',
        labels = output + 'predictive_labs_data_subset/CyTOF/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-{cell_num}_cells-{cell_selection}.tsv'
    output:
        prediction = output + 'predictive-labeling-benchmarking2/CyTOF/Init_{initial}-sel-{selection_procedure}-strat-{strat}-ALAlg-{AL_alg}-pred_alg-{pred_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-{neighbors}-resolution-{res}-seed-{s}-method-CyTOF-LDA-predictions-{cell_num}-cells-{cell_selection}.tsv'
    script:
        'cell-type-assignment/CyTOFLDA.R'
        
rule benchmark_predictive_labeling:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        active_learning_entropy = lambda wildcards: pred_lab_expand_normal_predictions_by_mod(wildcards.modality, "Active-Learning_entropy"),
        active_learning_maxp = lambda wildcards: pred_lab_expand_normal_predictions_by_mod(wildcards.modality, "Active-Learning_maxp"),
        pred_lab_files = lambda wildcards: get_pred_labeling(wildcards.modality),
        marker_seurat_files = lambda wildcards: pred_lab_expand_normal_predictions_by_mod(wildcards.modality, "MarkerSeurat-clustering"),
        NoMarker_seurat_files = lambda wildcards: pred_lab_expand_normal_predictions_by_mod(wildcards.modality, "NoMarkerSeurat-clustering"),
        random_files = lambda wildcards: pred_lab_expand_normal_predictions_by_mod(wildcards.modality, "random")
    output:
        acc = output + 'new/benchmark-predictive-labeling-{modality}-seed-{s}.tsv'
    log:
        output + 'logs/benchmark-predictive-labeling-{modality}-seed-{s}.log'
    resources:
        mem_mb=get_mem_mb
    script:
        'benchmarking/save-predictive-labeling-acc.R'

rule combine_predictive_labeling:
    input:
        accs = expand(output + 'new/benchmark-predictive-labeling-{{modality}}-seed-{s_val}.tsv', s_val = train_test_seeds)
    output:
        acc = output + 'new/pred2/benchmark-predictive-labeling-{modality}.tsv'
    script:
        'benchmarking/combine-pred-labeling-accs.R'

rule viz_benchmark_pred_labeling:
    input:
        acc = output + 'new/pred2/benchmark-predictive-labeling-{modality}.tsv'
    output:
        html = output + 'reports/pred2/benchmark-predictive-labeling-{modality}.html'
    params:
        output_dir = output + 'reports/pred2/'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/benchmark-predictive-labeling.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(acc = '{input.acc}'))\" "


rule identify_mislabelled:
    input:
        sce = 'data/{modality}/{modality}-train-seed-{seed}.rds'
    output:
        pred = output + 'identify_mislabelled/{modality}/Pred_alg-{pred_alg}-corr-10-seed-{seed}-cellNum-250-cells.tsv'
    resources:
        mem_mb=5000
    script:
        'predictive-labeling/detect-mislabelled.R'
