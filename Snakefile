
import pandas as pd
import numpy as np
import yaml
from itertools import chain

container: "al_eval_docker.sif" #"docker://mgeuenich/al_eval_docker"

configfile: 'config/config.yml'
output = 'output/' + config['version'] + '/'

AL_methods = ['multinom', 'rf']
selection_procedures = ['random', 'Seurat-clustering', 'Ranked-Active-Learning_entropy', 'Ranked-Active-Learning_maxp', 'Random-Active-Learning_entropy', 'Random-Active-Learning_maxp']
evaluation_methods_dict = {
    'scRNASeq': ['scmap-cluster', 'scmap-sc', 'singleR', 'Random-Forest', 'SVM-rejection', 'singleCellNet'], 
    'CyTOF': ['CyTOF-LDA', 'Random-Forest'],
    'snRNASeq': ['scmap-cluster', 'scmap-sc', 'singleR', 'Random-Forest', 'SVM-rejection', 'singleCellNet'],
    'scRNALung': ['scmap-cluster', 'scmap-sc', 'singleR', 'Random-Forest', 'SVM-rejection', 'singleCellNet'],
    'liverAtlas': ['scmap-cluster', 'scmap-sc', 'singleR', 'Random-Forest', 'SVM-rejection', 'singleCellNet'],
    'tabulaVasc': ['scmap-cluster', 'scmap-sc', 'singleR', 'Random-Forest', 'SVM-rejection', 'singleCellNet'],
    }
evaluation_methods = list(chain.from_iterable(list(evaluation_methods_dict.values())))
train_test_seeds = list(range(10))
annotators = ['GroundTruth']
modalities = ['scRNASeq', 'CyTOF', 'snRNASeq', 'scRNALung', 'liverAtlas', 'liverAtlas2', 'tabulaVasc']
data_splits = ['train', 'test']
initial_selections = ['ranking', 'random']

# Seurat params
Seurat_neighbors = [10,20,30]
Seurat_resolution = [0.4,0.8,1.2]

random_sets = ['set1']
cell_numbers = [100, 250, 500]
corruption_percentages = [0]
random_percentages = [0]


# Get markers for each cohort
with open(r'markers/scRNASeq.yml') as file:
    cell_types = yaml.full_load(file)

with open(r'markers/CyTOF.yml') as file:
    cytof_types = yaml.full_load(file)

with open(r'markers/snRNASeq.yml') as file:
    sn_cell_types = yaml.full_load(file)

with open(r'markers/scRNALung.yml') as file:
    sclung_cell_types = yaml.full_load(file)

with open(r'markers/liverAtlas.yml') as file:
    liver_cell_types = yaml.full_load(file)

with open(r'markers/tabulaVasc.yml') as file:
    vasc_cell_types = yaml.full_load(file)

# clean up marker names
scRNA_cell_types = list(cell_types['cell_types'].keys())
scRNA_cell_types_clean = []
for i in scRNA_cell_types:
    j = i.replace(' ', '-')
    scRNA_cell_types_clean.append(j)

cytof_cell_types = list(cytof_types['cell_types'].keys())
cytof_cell_types_clean = []
for i in cytof_cell_types:
    j = i.replace(' ', '-')
    cytof_cell_types_clean.append(j)

snRNASeq_cell_types = list(sn_cell_types['cell_types'].keys())
snRNASeq_cell_types_clean = []
for i in snRNASeq_cell_types:
    j = i.replace(' ', '-')
    snRNASeq_cell_types_clean.append(j)

scRNALung_cell_types = list(sclung_cell_types['cell_types'].keys())

liverAtlas_cell_types = list(liver_cell_types['cell_types'].keys())
liverAtlas_cell_types_clean = []
for i in liverAtlas_cell_types:
    j = i.replace(' ', '-')
    liverAtlas_cell_types_clean.append(j)

tabulaVasc_cell_types = list(vasc_cell_types['cell_types'].keys())
tabulaVasc_cell_types_clean = []
for i in tabulaVasc_cell_types:
    j = i.replace(' ', '-')
    tabulaVasc_cell_types_clean.append(j)

all_cell_types = {
    'scRNASeq': scRNA_cell_types_clean, 
    'CyTOF': cytof_cell_types_clean, 
    'snRNASeq': snRNASeq_cell_types_clean,
    'scRNALung': scRNALung_cell_types,
    'liverAtlas': liverAtlas_cell_types_clean,
    'tabulaVasc': tabulaVasc_cell_types_clean
}

# Needed for cell type removal in AL - cell type ID's need to match exactly
original_cell_types = {
    'scRNASeq': scRNA_cell_types, 
    'CyTOF': cytof_cell_types, 
    'snRNASeq': snRNASeq_cell_types,
    'scRNALung': scRNALung_cell_types,
    'liverAtlas': liverAtlas_cell_types,
    'tabulaVasc': tabulaVasc_cell_types
}

selection_expansion_dict = {
    'NoMarkerSeurat-clustering': {
        'initial': 'NA',
        'neighbors': Seurat_neighbors,
        'res': Seurat_resolution,
        'strategy': 'NA',
        'AL_alg': 'NA',
        'random_selection': [0],
        'corruption': [0]
    },
    'MarkerSeurat-clustering': {
        'initial': 'NA',
        'neighbors': Seurat_neighbors,
        'res': Seurat_resolution,
        'strategy': 'NA',
        'AL_alg': 'NA',
        'random_selection': [0],
        'corruption': [0]
    },
    'random': {
        'initial': 'NA',
        'neighbors': ['NA'],
        'res': ['NA'],
        'strategy': 'NA',
        'AL_alg': 'NA',
        'random_selection': [0],
        'corruption': [0]
    },
    'Active-Learning_entropy': {
        'initial': initial_selections,
        'neighbors': ['NA'],
        'res': ['NA'],
        'strategy': ['0.75_quant_entropy', '0.95_quant_entropy', 'highest_entropy'],
        'AL_alg': AL_methods,
        'random_selection': random_percentages,
        'corruption': corruption_percentages
    },
    'Active-Learning_maxp': {
        'initial': initial_selections,
        'neighbors': ['NA'],
        'res': ['NA'],
        'strategy': ['0.05_quant_maxp', '0.25_quant_maxp', 'lowest_maxp'],
        'AL_alg': AL_methods,
        'random_selection': random_percentages,
        'corruption': corruption_percentages
    }
}

#include: 'pipeline/process-data.smk'
#include: 'pipeline/cell-type-predictions.smk'
#include: 'pipeline/simulate-active-learning.smk'
include: 'pipeline/imbalance.smk'
include: 'pipeline/visualizations.smk'
# include: 'pipeline/imbalance2.smk'
#include: 'pipeline/rem-cell-type.smk'
# include: 'pipeline/predictive-labeling2.smk'
# include: 'pipeline/cell-type-similarity.smk'
# #include: 'pipeline/runtime.smk'
#include: 'pipeline/cell-type-sel-rem-markers.smk'
#include: 'pipeline/paper-figures.smk'

rule all:
    input:
        #process_data_output.values(),
        #[expand(x, s = train_test_seeds) for x in list(cell_type_predictions.values())],
        #active_learner.values(),
        viz.values(),
        # # imbalance.values(),
        # # imbalance2.values(),
        # rem_cell_type.values(),
        # pred_lab2.values(),
        # similarity.values(),
        #AR_remove_markers.values(),
        # final_figures.values()

