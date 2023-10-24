
final_figures = {
    'fig1': output + "paper-figures/figure1.pdf",
    # 'fig2': output + "paper-figures/lr-and-rf-heatmaps.pdf",
    # 'fig3_imb': output + "paper-figures/imbalance-main-figure.pdf",
    # 'fig3_acc': output + "paper-figures/overall-accuracies.pdf",
    'fig4': output + "paper-figures/rem-cell-type.pdf",
    # 'fig5': output + "paper-figures/pred-labelling.pdf",

    # # supplemental figures
    # 'ar_params': output + "paper-figures/Supp-AR-res.pdf"
}

# Figure 1
rule benchmarking_overview:
    input:
        cytof = "data/CyTOF/CyTOF-full.rds",
        scrna = "data/scRNASeq/scRNASeq-full.rds",
        snrna = "data/snRNASeq/snRNASeq-full.rds",
        scRNALung = "data/scRNALung/scRNALung-full.rds",
        liverAtlas = "data/liverAtlas/liverAtlas-full.rds",
        tabulaVasc = "data/tabulaVasc/tabulaVasc-full.rds",
        schematic = "illustrator-figures/benchmarking-schematic.ai"
    output:
        fig1 = output + "paper-figures/figure1.pdf"
    script:
        "figures/figure1.R"

# Figure 2
rule rf_vs_lr_and_ranking_vs_random:
    input:
        cytof = "data/CyTOF/CyTOF-train-seed-0.rds",
        cytof_AL_rand = expand(output + "data/CyTOF/Active-Learning_entropy/AL-batches-subset/Init_random-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        cytof_AL_rank = expand(output + "data/CyTOF/Active-Learning_entropy/AL-batches-subset/Init_ranking-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        scrna = "data/scRNASeq/scRNASeq-train-seed-0.rds",
        scrna_AL_rand = expand(output + "data/scRNASeq/Active-Learning_entropy/AL-batches-subset/Init_random-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        scrna_AL_rank = expand(output + "data/scRNASeq/Active-Learning_entropy/AL-batches-subset/Init_ranking-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        snrna = "data/snRNASeq/snRNASeq-train-seed-0.rds",
        snrna_AL_rand = expand(output + "data/scRNASeq/Active-Learning_entropy/AL-batches-subset/Init_random-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        snrna_AL_rank = expand(output + "data/scRNASeq/Active-Learning_entropy/AL-batches-subset/Init_ranking-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        snrna_lung = "data/scRNALung/scRNALung-train-seed-0.rds",
        snrna_lung_AL_rand = expand(output + "data/scRNALung/Active-Learning_entropy/AL-batches-subset/Init_random-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        snrna_lung_AL_rank = expand(output + "data/scRNALung/Active-Learning_entropy/AL-batches-subset/Init_ranking-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        liverAtlas = "data/liverAtlas/liverAtlas-train-seed-0.rds",
        liverAtlas_AL_rand = expand(output + "data/liverAtlas/Active-Learning_entropy/AL-batches-subset/Init_random-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        liverAtlas_AL_rank = expand(output + "data/liverAtlas/Active-Learning_entropy/AL-batches-subset/Init_ranking-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        tabulaVasc = "data/tabulaVasc/tabulaVasc-train-seed-0.rds",
        tabulaVasc_AL_rand = expand(output + "data/tabulaVasc/Active-Learning_entropy/AL-batches-subset/Init_random-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        tabulaVasc_AL_rank = expand(output + "data/tabulaVasc/Active-Learning_entropy/AL-batches-subset/Init_ranking-strat-highest_entropy-ALAlg-rf-rand_sel-0-corr-0-knn_neighbors-NA-resolution-NA-seed-{s}-100_cells.tsv", s = train_test_seeds),
        accs = expand(output + "results/overall-{mod}-benchmarking-accuracies.tsv", mod = modalities)
    output:
        rank_vs_rand_cell_num = output + "paper-figures/ranking-vs-random-cell-number.pdf",
        row_1_ranking_vs_random = output + "paper-figures/ranking-random-heatmap-row1.pdf",
        row_2_ranking_vs_random = output + "paper-figures/ranking-random-heatmap-row2.pdf",
        lr_rf_hm_1 = output + "paper-figures/lr-and-rf-heatmaps-row1.pdf",
        lr_rf_hm_2 = output + "paper-figures/lr-and-rf-heatmaps-row2.pdf",
        rf_lr_supp_1 = output + "paper-figures/rf-vs-lr-supplementary-1.pdf",
        rf_lr_supp_2 = output + "paper-figures/rf-vs-lr-supplementary-2.pdf",
        rand_rank_supp_1 = output + "paper-figures/random-vs-ranking-supplementary-1.pdf",
        rand_rank_supp_2 = output + "paper-figures/random-vs-ranking-supplementary-2.pdf"
    script:
        "figures/compare-lr-rf-random-and-ranking.R"

# Figure 3
rule imbalance_plot:
    input:
        cytof_acc = output + "imbalance/acc/imbalance-acc-CyTOF.tsv",
        scrna_acc = output + "imbalance/acc/imbalance-acc-scRNASeq.tsv",
        snrna_acc = output + "imbalance/acc/imbalance-acc-snRNASeq.tsv"
    output:
        main_fig = output + "paper-figures/imbalance-main-figure.pdf",
        sup_fig = output + "paper-figures/imbalance-supplementary.pdf"
    script:
        "figures/imbalance-plots.R"

rule overall_acc:
    input:
        accs = expand(output + "results/overall-{mod}-benchmarking-accuracies.tsv", mod = modalities),
    output:
        hm = output + "paper-figures/metric-correlation.pdf",
        overall_fig = output + "paper-figures/overall-accuracies.pdf"
    script:
        "figures/figure3.R"


# Figure 4
rule paper_rem_cell_types:
    input:
        cytof_summary = output + "figures/rem_cell_type_summary/entropy-CyTOF-summary.tsv",
        scrna_summary = output + "figures/rem_cell_type_summary/entropy-scRNASeq-summary.tsv",
        snrna_summary = output + "figures/rem_cell_type_summary/entropy-snRNASeq-summary.tsv",
        scrna_lung_summary = output + "figures/rem_cell_type_summary/entropy-scRNALung-summary.tsv",
        liver_summary = output + "figures/rem_cell_type_summary/entropy-liverAtlas-summary.tsv",
        vasc_summary = output + "figures/rem_cell_type_summary/entropy-tabulaVasc-summary.tsv",
        scrna = expand(output + "results/rem_cell_type/Init-sel-random-rem-{cell_line}-scRNASeq-highest_entropy-ALAlg-{al}-cell_num-20-seed-{s}.tsv", cell_line = ['JIMT1', 'AU565'], al = ['multinom', 'rf'], s = train_test_seeds),
        snrna = expand(output + "results/rem_cell_type/Init-sel-random-rem-{cell_type}-snRNASeq-highest_entropy-ALAlg-multinom-cell_num-20-seed-{s}.tsv", cell_type = ['Ductal', 'Endothelial', 'Schwann'], s = train_test_seeds),
        snrna_group = expand(output + "results/rem_cell_type_group/Init-sel-random-rem-schwann_l1-snRNASeq-highest_entropy-ALAlg-multinom-cell_num-20-seed-{s}.tsv", s = train_test_seeds),
        cytof = expand(output + "results/rem_cell_type/Init-sel-random-rem-{cell_type}-CyTOF-highest_entropy-ALAlg-multinom-cell_num-20-seed-{s}.tsv", cell_type = ['Classical Monocytes', 'CD8 T cells'], s = train_test_seeds),
        snrna_supp = expand(output + "results/rem_cell_type/Init-sel-random-rem-{cell_type}-snRNASeq-highest_entropy-ALAlg-rf-cell_num-20-seed-{s}.tsv", cell_type = ['Ductal', 'Endothelial', 'Schwann'], s = train_test_seeds),
        snrna_supp_group = expand(output + "results/rem_cell_type_group/Init-sel-random-rem-schwann_l1-snRNASeq-highest_entropy-ALAlg-rf-cell_num-20-seed-{s}.tsv", s = train_test_seeds),
        cytof_supp = expand(output + "results/rem_cell_type/Init-sel-random-rem-{cell_type}-CyTOF-highest_entropy-ALAlg-rf-cell_num-20-seed-{s}.tsv", cell_type = ['Classical Monocytes', 'CD8 T cells'], s = train_test_seeds),
    output:
        main = output + "paper-figures/rem-cell-type.pdf",
        supp = output + "paper-figures/supp-rem-cell-type.pdf"
    script:
        "figures/rem-cell-type.R"

# Figure 5
rule self_training:
    input:
        acc = output + "new/pred-labeling-accuracy.tsv",
        scrna = output + "new/pred2/benchmark-predictive-labeling-scRNASeq.tsv",
        snrna = output + "new/pred2/benchmark-predictive-labeling-snRNASeq.tsv",
        cytof = output + "new/pred2/benchmark-predictive-labeling-CyTOF.tsv",
        mislabelled_cytof = expand(output + "identify_mislabelled/CyTOF/Pred_alg-{al}-corr-10-seed-{s}-cellNum-250-cells.tsv", al = ['rf', 'multinom'], s = train_test_seeds),
        mislabelled_scrna = expand(output + "identify_mislabelled/scRNASeq/Pred_alg-{al}-corr-10-seed-{s}-cellNum-250-cells.tsv", al = ['rf', 'multinom'], s = train_test_seeds),
        mislabelled_snrna = expand(output + "identify_mislabelled/snRNASeq/Pred_alg-{al}-corr-10-seed-{s}-cellNum-250-cells.tsv", al = ['rf', 'multinom'], s = train_test_seeds),
    output:
        supp = output + "paper-figures/Supp-pred-labelling-acc.pdf",
        main = output + "paper-figures/pred-labelling.pdf",
        sup_cytof = output + "paper-figures/Supp-CyTOF-f1-improvement.pdf",
        sup_scrna = output + "paper-figures/Supp-scRNASeq-f1-improvement.pdf",
        sup_snrna = output + "paper-figures/Supp-snRNASeq-f1-improvement.pdf"
    script:
        "figures/pred-labeling.R"


#### Supplemental figs
# ar params
rule AR_params:
    input:
        cytof_acc = output + "results/overall-CyTOF-benchmarking-accuracies.tsv",
        scrna_acc = output + "results/overall-scRNASeq-benchmarking-accuracies.tsv",
        snrna_acc = output + "results/overall-snRNASeq-benchmarking-accuracies.tsv",
        scrna_lung = output + "results/overall-scRNALung-benchmarking-accuracies.tsv",
        liver = output + "results/overall-liverAtlas-benchmarking-accuracies.tsv",
        vasc = output + "results/overall-tabulaVasc-benchmarking-accuracies.tsv"
    output:
        res = output + "paper-figures/Supp-AR-res.pdf",
        knn = output + "paper-figures/Supp-AR-knn.pdf"
    script:
        "figures/AR-params.R"