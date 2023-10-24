packages <- c(
    'tidyverse',
    'RColorBrewer',
    'patchwork',
    'caret',
    'SingleCellExperiment',
    'splitstackshape',
    'fabricatr',
    'scmap',
    'SingleR',
    'Seurat',
    'yaml',
    'glue',
    'scater',
    'yardstick',
    'data.table',
    'scales',
    'ComplexHeatmap',
    'circlize',
    'viridis',
    'ggalluvial',
    'magick',
    'ggpubr',
    'broom'
)

ip <- installed.packages()
packages <- packages[!(packages %in% rownames(ip))]

BiocManager::install(packages)