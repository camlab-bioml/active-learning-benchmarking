
rem_cell_type_reports = []
rctreport = [expand(output + 'results/rem_cell_type/Init-sel-{initial}-rem-{rem_cell_type}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}-seed-{s}.tsv',
    initial = initial_selections,
    rem_cell_type = original_cell_types[modality],
    modality = [modality],
    AL_alg = AL_methods,
    strat = ['highest_entropy', 'lowest_maxp'],
    num = [20, 50],
    s = train_test_seeds)
    for modality in original_cell_types.keys()]

for element in rctreport:
    rem_cell_type_reports.extend(element)

group_rem_dict = {
    'CyTOF': ['B_cells', 'monocytes', 'T', 'T_NK'],
    'snRNASeq': ['tumour', 'schwann_l1', 'schwann_l2'],
    'CyTOFExpanded': ['CD8', 'pDCs', 'B_cells_all', 'B_cells_IgD_neg', 'B_and_prog']
}

rem_cell_type = {
    'plot': expand(output + 'figures/rem_cell_type/Selected_cells-Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf',
        initial = initial_selections, modality = modalities, strat = ['highest_entropy', 'lowest_maxp'], AL_alg = AL_methods, num = [20, 50]),
    'probs': expand(output + 'figures/rem_cell_type_prob/Selected_cells-Init-sel-{initial}-{modality}-ALAlg-{AL_alg}-cell_num-{num}-rem_celltype-{rem_cell_type}-seed-{s}.pdf',
        initial = ['random'], modality = ['scRNASeq'], AL_alg = AL_methods, num = [20], rem_cell_type = original_cell_types['scRNASeq'], s = [0]),
    'probs_sn': expand(output + 'figures/rem_cell_type_prob/Selected_cells-Init-sel-{initial}-{modality}-ALAlg-{AL_alg}-cell_num-{num}-rem_celltype-{rem_cell_type}-seed-{s}.pdf',
        initial = ['random'], modality = ['snRNASeq'], AL_alg = AL_methods, num = [20], rem_cell_type = original_cell_types['snRNASeq'], s = [0]),

    'rem_group': expand(output + 'figures/rem_cell_type_group/Selected_cells-Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf',
        initial = initial_selections, modality = ['CyTOF', 'snRNASeq'], strat = ['highest_entropy', 'lowest_maxp'], AL_alg = AL_methods, num = [20]),
    
    'debug_rf': expand(output + 'figures/rem_cell_type_debug/entropy-rem_cell_type-{rem}.pdf', rem = original_cell_types['snRNASeq']),
    'summary': expand(output + 'figures/rem_cell_type_summary/uncertainty-summary-{modality}.png', modality = modalities)
}


def get_mem_mb(wildcards, attempt):
    return attempt * 3000 + 2000

rule rem_cell_type:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    resources:
        mem_mb=get_mem_mb
    output:
        tsv = output + 'results/rem_cell_type/Init-sel-{initial}-rem-{rem_cell_type}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}-seed-{s}.tsv'
    script:
        "rem-cell-type-from-training/rem-cell-type.R"

def get_entropy_values(mod, cell_type_markers = original_cell_types):
    base_files = output + 'results/rem_cell_type/Init-sel-{{initial}}-rem-{rem_cell_type}-{modality}-{{strat}}-ALAlg-{{AL_alg}}-cell_num-20-seed-{{s}}.tsv'
    files  = expand(base_files, rem_cell_type = cell_type_markers[mod], modality = mod)

    rf_maxp_files = expand(files, initial = ['ranking', 'random'], strat = ['lowest_maxp'], AL_alg = ['rf'], s = train_test_seeds)
    rf_entr_files = expand(files, initial = ['ranking', 'random'], strat = ['highest_entropy'], AL_alg = ['rf'], s = train_test_seeds)
    lr_maxp_files = expand(files, initial = ['ranking', 'random'], strat = ['lowest_maxp'], AL_alg = ['multinom'], s = train_test_seeds)
    lr_entr_files = expand(files, initial = ['ranking', 'random'], strat = ['highest_entropy'], AL_alg = ['multinom'], s = train_test_seeds)
    
    return rf_maxp_files + rf_entr_files + lr_maxp_files + lr_entr_files

rule summarize_entropies:
    input:
        entropy_files = lambda wildcards: get_entropy_values(wildcards.modality)
    output:
        entropies = output + 'figures/rem_cell_type_summary/uncertainty-summary-{modality}.png',
        example_plot = output + 'figures/rem_cell_type_summary/entropy-{modality}-summary.tsv'
    script:
        'rem-cell-type-from-training/summarize-entropies.R'


def get_individual_cell_types(mod, ct_group):
    if mod == "CyTOF":
        if ct_group == "B_cells":
            cts = ['IgD- IgMpos B cells', 'IgM- IgD- B-cells', 'IgDpos IgMpos B cells']
        elif ct_group == "monocytes":
            cts = ['Intermediate Monocytes', 'Eosinophils', 'Classical Monocytes']
        elif ct_group == "T":
            cts = ['CD8 T cells', "CD4 T cells"]
        elif ct_group == "T_NK":
            cts = ['CD8 T cells', "CD4 T cells", 'NKT cells', 'NK cells']

    if mod == "CyTOFExpanded":
        if ct_group == "CD8":
            cts = ['CD8 T cells']
        elif ct_group == "pDCs":
            cts = ['pDCs']
        elif ct_group == "B_cells_all":
            cts = ['IgD- IgMpos B cells', 'IgM- IgD- B-cells', 'B-cell Frac A-C (pro-B cells)']
        elif ct_group == "B_cells_IgD_neg":
            cts = ['IgD- IgMpos B cells', 'IgM- IgD- B-cells']
        elif ct_group == "B_and_prog":
            cts = ['IgD- IgMpos B cells', 'IgM- IgD- B-cells', 'B-cell Frac A-C (pro-B cells)', 'CMP', 'GMP']
    
    
    if mod == "scRNASeq":
        if ct_group == "T_NK":
            cts = ['CD4+ T cell', 'Cytotoxic T cell', 'Natural killer cell']
    
    if mod == "snRNASeq":
        if ct_group == "tumour":
            cts = ['Atypical_Ductal', 'Tumor']
        elif ct_group == "schwann_l1":
            cts = ['Endothelial', 'Schwann']
        elif ct_group == "schwann_l2":
            cts = ['Endothelial', 'Schwann', 'SmoothMuscle', 'Fibroblast']
    
    return cts

rule rem_group_cell_types:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    params:
        rem_cell_type_list = lambda wildcards: get_individual_cell_types(wildcards.modality, wildcards.rem_cell_type_group)
    resources:
        mem_mb=10000
    output:
        tsv = output + 'results/rem_cell_type_group/Init-sel-{initial}-rem-{rem_cell_type_group}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}-seed-{s}.tsv'
    script:
        "rem-cell-type-from-training/rem-cell-type.R"

rule combine_uncertainties:
    input:
        tsv = rem_cell_type_reports
    output:
        tsv = output + 'results/rem_cell_type/Rem-cell-types-combined-uncertanties.tsv'
    script:
        'rem-cell-type-from-training/combine-uncertainties.R'

def get_uncertainties_by_modality(mod, mod_cts = original_cell_types, seeds = train_test_seeds):
    cts = mod_cts[mod]
    f = expand(output + 'results/rem_cell_type/Init-sel-{{initial}}-rem-{rem_cell_type}-{{modality}}-{{strat}}-ALAlg-{{AL_alg}}-cell_num-{{num}}-seed-{s}.tsv',
        rem_cell_type = cts, s = seeds)
    
    return f


rule plot_proportion_selected:
    input:
        tsvs = lambda wildcards: get_uncertainties_by_modality(wildcards.modality)
    output:
        plot = output + 'figures/rem_cell_type/Selected_cells-Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf'
    script:
        'rem-cell-type-from-training/plot-proportion-selected.R'


def get_uncertainties_by_modality_group(mod):
    if mod == "CyTOFExpanded":
        f = expand(output + 'results/rem_cell_type_group/Init-sel-{{initial}}-rem-{rem_cell_type_group}-{{modality}}-{{strat}}-ALAlg-{{AL_alg}}-cell_num-{{num}}-seed-0.tsv',
            rem_cell_type_group = group_rem_dict[mod])
    else:
        f = expand(output + 'results/rem_cell_type_group/Init-sel-{{initial}}-rem-{rem_cell_type_group}-{{modality}}-{{strat}}-ALAlg-{{AL_alg}}-cell_num-{{num}}-seed-{s}.tsv',
            rem_cell_type_group = group_rem_dict[mod], s = train_test_seeds)
    return f


rule entropy_maxp_group_rem:
    input:
        tsvs = lambda wildcards: get_uncertainties_by_modality_group(wildcards.modality)
    output:
        plot = output + 'figures/rem_cell_type_group/Selected_cells-Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf'
    script:
        'rem-cell-type-from-training/plot-proportion-selected.R'

rule understand_probabilities:
    input:
        markers = 'markers/{modality}.yml',
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    resources: mem_mb=5000
    output:
        pdf = output + 'figures/rem_cell_type_prob/Selected_cells-Init-sel-{initial}-{modality}-ALAlg-{AL_alg}-cell_num-{num}-rem_celltype-{rem_cell_type}-seed-{s}.pdf'
    script:
        'rem-cell-type-from-training/understand-probs.R'


rule debug_rf:
    input:
        sce = "data/snRNASeq/snRNASeq-train-seed-0.rds",
        test = "data/snRNASeq/snRNASeq-test-seed-0.rds"
    output:
        probs_entropy = output + 'figures/rem_cell_type_debug/entropy-rem_cell_type-{rem}.pdf'
    script:
        'rem-cell-type-from-training/debug-random-forest.R'

# rule find_doublets:
#     input:
#         sce = 'data/scRNASeq/scRNASeq-full.rds'
#     params:
#         output_dir = output + 'reports/'
#     output:
#         html = output + 'reports/scRNASeq-doublet-id-{sample}.html',
#         doublets = output + 'results/doublet-id/scRNASeq-doublet-id-{sample}.tsv'
#     shell:
#         "Rscript -e \"Sys.setenv(RSTUDIO_PANDOC='/u/ext.mgeuenich/pandoc/pandoc-2.11.4/bin/');"
#         "rmarkdown::render('pipeline/rem-cell-type-from-training/DoubletFinder.Rmd', output_dir = '{params.output_dir}', output_file = '{output.html}', "
#         "params = list(sce = '{input.sce}', method = '{wildcards.sample}', res = '{output.doublets}'))  \" "

# rule plot_doublet_entropy:
#     input:
#         doublets1 = output + 'results/doublet-id/scRNASeq-doublet-id-10x Chromium (v2) A.tsv',
#         doublets2 = output + 'results/doublet-id/scRNASeq-doublet-id-10x Chromium (v2) B.tsv',
#         entropies = lambda wildcards: get_uncertainties_by_modality(wildcards.modality)
#     output:
#         pdf = output + 'figures/doublet-rem-cell-type/Init-sel-{initial}-{modality}-{strat}-ALAlg-{AL_alg}-cell_num-{num}.pdf'
#     script:
#         'rem-cell-type-from-training/doublet-finder-entropies.R'
