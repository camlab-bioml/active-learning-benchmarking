def create_AL_iterations_assignments(mod):
    if(mod == 'CyTOF'):
        methods = evaluation_methods_dict['CyTOF']
    else:
        methods = evaluation_methods_dict['scRNASeq']

    return(expand(output + 'rare-subtype-benchmarking/{{modality}}-Active-Learning-annotator-GroundTruth-max_dim-NA-resolution-NA-iterations_set-{set}-{method}-predictions.tsv',
           set = list(range(training_set_AL_simulation)), method = methods))

def get_predictions(mod, seed):
    if mod == 'scRNASeq':
        f = expand(cell_type_predictions['scRNASeq'], s = seed)
    elif mod == 'snRNASeq':
        f = expand(cell_type_predictions['snRNASeq'], s = seed)
    elif mod == 'scRNALung':
        f = expand(cell_type_predictions['scRNALung'], s = seed)
    elif mod == 'CyTOF':
        f = expand(cell_type_predictions['CyTOF'], s = seed)
    elif mod == 'liverAtlas':
        f = expand(cell_type_predictions['liverAtlas'], s = seed)
    elif mod == 'liverAtlas2':
        f = expand(cell_type_predictions['liverAtlas2'], s = seed)
    elif mod == 'tabulaVasc':
        f = expand(cell_type_predictions['tabulaVasc'], s = seed)
    
    return (f)

alluvials = []
alluv = [expand(output + 'figures/gt_predictions/{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}_rand_sel-0-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf',
            modality = ['liverAtlas', 'scRNASeq'],
            selection_procedure = [select], 
            initial = small_selection_strat_dict[select]['initial'],
            strat = small_selection_strat_dict[select]['strategy'],
            AL_alg = small_selection_strat_dict[select]['AL_alg'],
            #rand = small_selection_strat_dict[select]['random_selection'],
            corruption = small_selection_strat_dict[select]['corruption'],
            knn = small_selection_strat_dict[select]['neighbors'], 
            res = small_selection_strat_dict[select]['res'], 
            method = evaluation_methods_dict['scRNASeq'],
            cell_num = [100, 250, 500]) 
            for select in small_selection_strat_dict.keys()]
for element in alluv:
    alluvials.extend(element)

viz = {
    # 'benchmark': expand(output + 'reports/overall-{modality}-{AL_alg}-benchmarking-cells.html', modality = modalities, AL_alg = AL_methods),
    #'cross_cohort': output + 'reports/cross-cohort-eval.html',
    'f1_AL_performance': expand(output + 'figures/AL_f1/f1-performance-{plot_num}.pdf', plot_num = [1, 2]),
    #'alluvs': alluvials
}

rule plot_AL_accuracy:
    input:
        maxp_accs = expand(output + 'results/AL_f1/{modality}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-F1-score.tsv',
            modality = modalities, initial = selection_expansion_dict['Active-Learning_maxp']['initial'], strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], 
            AL_alg = ['rf', 'multinom'], rand = [0], corrupt = [0, 1], s = train_test_seeds),
        entr_accs = expand(output + 'results/AL_f1/{modality}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-F1-score.tsv',
            modality = modalities, initial = selection_expansion_dict['Active-Learning_entropy']['initial'], strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], 
            AL_alg = ['rf', 'multinom'], rand = [0], corrupt = [0, 1], s = train_test_seeds)
    output:
        pdf1 = output + 'figures/AL_f1/f1-performance-1.pdf',
        pdf2 = output + 'figures/AL_f1/f1-performance-2.pdf'
    script:
        'visualize/visualize-AL-classifier-accuracy.R'

rule cross_cohort_overview:
    input:
        accs = expand(output + 'results/overall-{modality}-benchmarking-accuracies.tsv', modality = modalities)
    output:
        html = output + 'reports/cross-cohort-eval.html'
    params:
        output_dir = output + 'reports/'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/Cross-cohort-eval.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(accs = '{input.accs}'))\" "

rule visualize_selected_cells:
    input:
        sel_cells = expand(output + 'data/{{modality}}/{{selection}}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{cell_num}_cells.tsv', cell_num = cell_numbers, s = train_test_seeds)
    output:
        barplot = output + 'figures/selected_cells/{modality}/{selection}/Init_NA-strat-NA-ALAlg-NA-rand_sel-0-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-selected-cells.R'

rule visualize_selected_cells_AL_corrupt:
    input:
        sel_cells = expand(output + 'data/{{modality}}/{{selection}}/AL-batches-subset/Init_{initial}-strat-{{strat}}-ALAlg-{{al}}-rand_sel-0-corr-{cor}-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{cell_num}_cells.tsv', 
            initial = initial_selections, cell_num = cell_numbers, s = train_test_seeds, cor = corruption_percentages)
    params:
        cor_or_rand = "cor"
    output:
        barplot = output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-0-corr-all-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-selected-cells.R'

rule visualize_selected_cells_AL_rand:
    input:
        sel_cells = expand(output + 'data/{{modality}}/{{selection}}/AL-batches-subset/Init_{initial}-strat-{{strat}}-ALAlg-{{al}}-rand_sel-{rand}-corr-0-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{cell_num}_cells.tsv', 
            initial = initial_selections, cell_num = cell_numbers, s = train_test_seeds, rand = random_percentages)
    params:
        cor_or_rand = "rand"
    output:
        barplot = output + 'figures/selected_cells/{modality}/{selection}/Init_all-strat-{strat}-ALAlg-{al}-rand_sel-all-corr-0-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-selected-cells.R'

def get_gt_predicted_input_files(mod, pred_labeling, imbalance = False):
    if pred_labeling:
        f = expand(output + 'predictive-labeling-benchmarking/{{modality}}/Init_{{initial}}-sel-{{selection_procedure}}-strat-{{strat}}-ALAlg-{{AL_alg}}-pred_alg-{{pred_lab_alg}}-rand_sel-{{rand}}-corr-{{corruption}}-knn_neighbors-{{knn}}-resolution-{{res}}-seed-{s}-{method}-predictions-{cell_num}-cells-{{pred_lab_selection}}.tsv', 
                method = evaluation_methods_dict[mod],
                s = train_test_seeds,
                cell_num = cell_numbers)
    elif imbalance:
        f = expand(output + 'imbalance/rare-subtype-benchmarking/{{bal}}-{{similarity}}/Init-{{initial}}-{{modality}}-sel-{{selection_procedure}}-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corruption}}-knn_neighbors-{{knn}}-resolution-{{res}}-{method}-predictions-seed-{s}-{cell_num}-cells.tsv', 
                method = evaluation_methods_dict[mod],
                s = train_test_seeds,
                cell_num = 100)
    else:
        f = expand(output + 'rare-subtype-benchmarking/Init_{{initial}}-{{modality}}-sel_{{selection_procedure}}-strat-{{strat}}-ALAlg-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corruption}}-knn_neighbors-{{knn}}-resolution-{{res}}-{method}-predictions-seed-{s}-{cell_num}-cells.tsv',
                method = evaluation_methods_dict[mod],
                s = train_test_seeds,
                cell_num = cell_numbers)
    return f

rule compare_gt_to_predicted:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        predicted = lambda wildcards: get_gt_predicted_input_files(wildcards.modality, False)
    output:
        pdf = output + 'figures/gt_predictions/{modality}/sel_{selection_procedure}/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}_rand_sel-{rand}-corr-{corruption}-knn_neighbors-{knn}-resolution-{res}.pdf'
    script:
        'visualize/plot-gt-against-predicted.R'


rule overall_benchmark:
    input:
        sce = 'data/{modality}/{modality}-full.rds',
        predictions = lambda wildcards: get_predictions(wildcards.modality, wildcards.s)
    resources:
        mem_mb=50000
    log:
        output + 'logs/benchmark-predictive-labeling-{modality}-seed-{s}.log'
    output:
        acc = output + 'results/overall-{modality}-benchmarking-accuracies-seed-{s}.tsv'
    script:
        'benchmarking/save-acc-overall-benchmarking.R'

rule combine_overall_benchmark:
    input:
        tsvs = expand(output + 'results/overall-{{modality}}-benchmarking-accuracies-seed-{s}.tsv', s = train_test_seeds)
    output:
        acc = output + 'results/overall-{modality}-benchmarking-accuracies.tsv'
    script:
        'benchmarking/combine-acc-overall-benchmarking.R'

rule visualize_overall_benchmark:
    input:
        acc = output + 'results/overall-{modality}-benchmarking-accuracies.tsv'
    output:
        html = output + 'reports/overall-{modality}-{AL_alg}-benchmarking-cells.html'
    params:
        output_dir = output + 'reports/'
    shell:
        "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/home/campbell/share/software/pandoc-2.9.2.1/bin');"
        "rmarkdown::render('pipeline/benchmarking/benchmarking.Rmd', output_file = '{output.html}', output_dir = '{params.output_dir}', "
        "params = list(acc = '{input.acc}', al_alg = '{wildcards.AL_alg}'))\" "
