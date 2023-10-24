
modalities = ['scRNASeq', 'CyTOF', 'snRNASeq', 'scRNALung', 'tabulaLiver', 'tabulaVasc']
train_test_seeds = list(range(10))
output = 'output/v8/'
runtime = {
    #'runtime_vals': expand(output + 'results/runtime/timing-{modality}-seed-{s}.csv', modality = modalities, s = train_test_seeds),
}


rule runtime:
    input:
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds',
        markers = 'markers/{modality}.yml',
    output:
        timing = output + 'results/runtime/timing-{modality}-seed-{s}.csv',
    script:
        'runtime/runtime.R'

rule combine_runtime:
    input:
        expand(output + 'results/runtime/timing-{modality}-seed-{s}.csv', modality = modalities, s = train_test_seeds),
    output:
        'all.txt'
    shell:
        'test'