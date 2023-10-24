

similarity = {
    'heatmaps': expand(output + 'figures/cell-type-similarity/similarity-heatmap-{modality}-seed-{s}.pdf', modality = modalities, s = train_test_seeds),
    #'cellline': output + 'figures/cell-type-similarity/similarity-heatmap-scRNASeq-cellLines.pdf',
    #'tsvs': expand(output + 'results/cell-type-similarity/similarity-{modality}-seed-{s}.tsv', modality = modalities, s = train_test_seeds),
    #'mean': expand(output + 'results/cell-type-similarity/average-similarity-{modality}.tsv', modality = modalities)
}

rule calculate_cosine_dist:
    input:
        sce = 'data/{modality}/{modality}-train-seed-{s}.rds'
    output:
        heatmap = output + 'figures/cell-type-similarity/similarity-heatmap-{modality}-seed-{s}.pdf',
        tsv = output + 'results/cell-type-similarity/similarity-{modality}-seed-{s}.tsv'
    script:
        'cell-type-similarity/calculate-cosine-distance.R'

rule cell_line_dists:
    input:
        sce = 'data/scRNA-cell-lines/scRNASeq-cellLines.rds'
    params:
        pca = 'scater' 
    resources: mem_mb=100000
    output:
        heatmap = output + 'figures/cell-type-similarity/similarity-heatmap-scRNASeq-cellLines.pdf',
        tsv = output + 'results/cell-type-similarity/similarity-scRNASeq-cellLines.tsv'
    script:
        'cell-type-similarity/calculate-cosine-distance.R'

rule summarize_cosine_dist_across_seeds:
    input:
        tsvs = expand(output + 'results/cell-type-similarity/similarity-{modality}-seed-{s}.tsv',
            modality = modalities, s = train_test_seeds)
    resources: mem_mb=10000
    output:
        avg_dist_sc = output + 'results/cell-type-similarity/average-similarity-scRNASeq.tsv',
        avg_dist_sn = output + 'results/cell-type-similarity/average-similarity-snRNASeq.tsv',
        avg_dist_cy = output + 'results/cell-type-similarity/average-similarity-CyTOF.tsv',
        avg_dist_scLung = output + 'results/cell-type-similarity/average-similarity-scRNALung.tsv'
    script:
        'cell-type-similarity/summarize-cosine-distances.R'