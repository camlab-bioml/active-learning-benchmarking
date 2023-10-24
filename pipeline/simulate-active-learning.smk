def get_entropy_files(AL_type):
    if AL_type == "Active-Learning_entropy":
        f = expand(output + 'data/{{modality}}/uncertainties/Active-Learning_entropy/Init-{{initial}}-strat-{strat}-al-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-seed-{{s}}.tsv',
            strat = selection_expansion_dict['Active-Learning_entropy']['strategy'])
    elif AL_type == 'Active-Learning_maxp':
        f = expand(output + 'data/{{modality}}/uncertainties/Active-Learning_maxp/Init-{{initial}}-strat-{strat}-al-{{AL_alg}}-rand_sel-{{rand}}-corr-{{corrupt}}-seed-{{s}}.tsv',
            strat = selection_expansion_dict2['Active-Learning_maxp']['strategy'])
    return f

active_learner = {
    'maxp_box_rand_0': expand(output + 'figures/uncertainty-{modality}-Init-{initial}-strat-{AL_type}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.pdf',
        modality = modalities, AL_type = ['Active-Learning_maxp', 'Active-Learning_entropy'], initial = initial_selections, 
        AL_alg = AL_methods,
        rand = 0, 
        corrupt = corruption_percentages,
        s = train_test_seeds),
    'maxp_box_cor_0': expand(output + 'figures/uncertainty-{modality}-Init-{initial}-strat-{AL_type}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.pdf',
        modality = modalities, AL_type = ['Active-Learning_maxp', 'Active-Learning_entropy'], initial = initial_selections, 
        AL_alg = AL_methods,
        rand = random_percentages, 
        corrupt = 0,
        s = train_test_seeds),

    'subset_al_entropy_rand_0': expand(output + 'data/{modality}/{AL_type}/AL-batches-subset/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-NA-resolution-NA-seed-{s}-{subset_val}_cells.tsv', 
        modality = modalities, initial = selection_expansion_dict['Active-Learning_entropy']['initial'], AL_type = ['Active-Learning_entropy'], strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], 
        AL_alg = selection_expansion_dict['Active-Learning_entropy']['AL_alg'],
        rand = 0, corrupt = selection_expansion_dict['Active-Learning_entropy']['corruption'], 
        subset_val = cell_numbers, s = train_test_seeds),

    'subset_al_entropy_corr_0': expand(output + 'data/{modality}/{AL_type}/AL-batches-subset/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-NA-resolution-NA-seed-{s}-{subset_val}_cells.tsv', 
        modality = modalities, initial = selection_expansion_dict['Active-Learning_entropy']['initial'], AL_type = ['Active-Learning_entropy'], strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], 
        AL_alg = selection_expansion_dict['Active-Learning_entropy']['AL_alg'],
        rand = selection_expansion_dict['Active-Learning_entropy']['random_selection'], corrupt = 0, 
        subset_val = cell_numbers, s = train_test_seeds),

    'f1_scores_entropy': expand(output + 'results/AL_f1/{modality}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-F1-score.tsv',
        modality = modalities, initial = selection_expansion_dict['Active-Learning_entropy']['initial'], strat = selection_expansion_dict['Active-Learning_entropy']['strategy'], AL_alg = ['rf', 'multinom'], rand = [0], corrupt = [0, 1], s = train_test_seeds),
    'f1_scores_maxp': expand(output + 'results/AL_f1/{modality}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-F1-score.tsv',
        modality = modalities, initial = selection_expansion_dict['Active-Learning_maxp']['initial'], strat = selection_expansion_dict['Active-Learning_maxp']['strategy'], AL_alg = ['rf', 'multinom'], rand = [0], corrupt = [0, 1], s = train_test_seeds)
}

rule Create_active_learning_ground_truth:
    input:
        markers = 'markers/{modality}.yml',
        expression = 'data/{modality}/{modality}-train-seed-{s}.rds'
    params:
        max_cell_num = max(cell_numbers)
    priority: 10
    resources:
        mem_mb=10000
    output:
        assignments = output + 'data/{modality}/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv',
        entropy = output + 'data/{modality}/uncertainties/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv'
    script:
        'cell-type-assignment/simulate-active-learner.R'

selection_expansion_dict2 = {
    'Active-Learning_maxp': {
       'initial': initial_selections,
       'neighbors': ['NA'],
               'res': ['NA'],
               'strategy': ['0.05_quant_maxp', '0.25_quant_maxp', 'lowest_maxp'],
               'initial_sel': initial_selections,
               'AL_alg': AL_methods,
               'random_selection': random_percentages,
               'corruption': corruption_percentages}
}

rule visualize_entropies:
    input:
        f = lambda wildcards: get_entropy_files(wildcards.AL_type)
    output:
        box = output + 'figures/uncertainty-{modality}-Init-{initial}-strat-{AL_type}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.pdf'
    script:
        'visualize/plot-entropies.R'

rule create_AL_training_batches:
    input:
        assignment = output + 'data/{modality}/{AL_type}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}.tsv'
    output:
        split = output + 'data/{modality}/{AL_type}/AL-batches-subset/Init_{initial}-strat-{strat}-ALAlg-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-knn_neighbors-NA-resolution-NA-seed-{s}-{subset_val}_cells.tsv'
    resources: smallfile=1
    script:
        'cell-type-assignment/subset-simulated-active-learner-n-cells.R'

rule check_AL_f1_score:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    resources:
        mem_mb=8000
    output:
        tsv = output + 'results/AL_f1/{modality}/Init-{initial}-strat-{strat}-al-{AL_alg}-rand_sel-{rand}-corr-{corrupt}-seed-{s}-F1-score.tsv'
    script:
        'cell-type-assignment/active-learning-accuracy.R'
