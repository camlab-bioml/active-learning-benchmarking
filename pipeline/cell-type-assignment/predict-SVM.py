import numpy as np
import pandas as pd
import yaml
from itertools import chain
import joblib


## Load everything 
model = joblib.load(open(snakemake.input['model'], 'rb'))
expression = pd.read_csv(snakemake.input['test'], header = 0, sep = '\t')

m = list(model.feature_names_in_) + ['cell_type', 'cell_id']

expression = expression[m]

true_labels = expression['cell_type']
cell_ids = expression['cell_id']
expression = expression.drop(['cell_type', 'cell_id'], axis = 1)

# predict probabilities and labels everything under 0.7 unassigned
predictions = model.predict_proba(expression)
predictions = pd.DataFrame(predictions, columns=model.classes_)

# create column with highest probability and other column with corresponding class
predictions['max_prob'] = predictions.max(axis=1)
predictions['class'] = predictions.idxmax(axis=1)

# overwrite class with 'unassigned' if max_prob is below 0.7
predictions.loc[predictions['max_prob'] < 0.7, 'class'] = 'unassigned'

predicted = predictions['class'].to_list()

if 'cell_num' in dict(snakemake.wildcards).keys():
    cell_num = snakemake.wildcards['cell_num']
else:
    cell_num = 'NA'

output = pd.DataFrame({'cell_id': cell_ids,
                       'predicted_cell_type': predicted,
                       'prediction_params': 'SVM-rejection-knn-' + snakemake.wildcards['neighbors'] + '-res-' + snakemake.wildcards['res'] + '-cell_numbers-' + cell_num + '-randomSelection-' + snakemake.wildcards['rand'] + '-corrupted-' + snakemake.wildcards['corrupt'] + '-Init-' + snakemake.wildcards['initial'] + '-seed-' + snakemake.wildcards['s'],
                       'selection_procedure': snakemake.wildcards['selection_procedure'] + '-strategy-' + snakemake.wildcards['strat'] + '-ALAlg-' + snakemake.wildcards['AL_alg'],
                       'modality': snakemake.wildcards['modality']})

if 'cell_selection' in dict(snakemake.wildcards).keys():
    output['cell_selection'] = snakemake.wildcards['cell_selection']
else:
    output['cell_selection'] = 'NA'

if 'similarity' in dict(snakemake.wildcards).keys():
    output['similarity'] = snakemake.wildcards['bal'] + '-' + snakemake.wildcards['similarity']

if 'similarity' not in dict(snakemake.wildcards).keys() and 'bal' in dict(snakemake.wildcards).keys():
    output['balance'] = snakemake.wildcards['bal']

if 'rem_percentage' in dict(snakemake.wildcards).keys():
    output['rem_percentage'] = snakemake.wildcards['rem_percentage']

if 'cell_selection' in dict(snakemake.wildcards).keys():
    output['pred_cells'] = snakemake.wildcards['cell_selection']

output.to_csv(snakemake.output['predictions'], sep = '\t', index = False)