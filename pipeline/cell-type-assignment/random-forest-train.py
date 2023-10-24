import numpy as np
import pandas as pd
import random

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
import joblib


# ## Data processing
expression = pd.read_csv(snakemake.input['train'], sep = "\t", header = 0)
annotation = pd.read_csv(str(snakemake.input['annotation']), sep = '\t', header = 0)

if 'cell_selection' in dict(snakemake.wildcards).keys():
    expression = pd.merge(expression, annotation, on = 'cell_id')
    X_train = expression.drop(['cell_type_y', 'entropy', 'labeling', 'cell_type_x', 'cell_id', 'corrupted_cell_type', 'iteration', 'method', 'cell_num', 'params'], axis = 1)
    y_train = expression['cell_type_y']
else:
    expression = pd.merge(expression, annotation, on = 'cell_id')

    if snakemake.wildcards['selection_procedure'] == 'Active-Learning_entropy' or snakemake.wildcards['selection_procedure'] == 'Active-Learning_maxp':
        expression = expression.drop(['iteration', 'cell_num', 'corrupted_cell_type'], axis = 1)

    if snakemake.wildcards['selection_procedure'] == 'random':
        expression = expression.drop(['set','params'], axis = 1)

    if snakemake.wildcards['selection_procedure'] == 'NoMarkerSeurat-clustering' or snakemake.wildcards['selection_procedure'] == 'MarkerSeurat-clustering':
        expression = expression.drop(['params'], axis = 1)

    X_train = expression.drop(['cell_type_y', 'method', 'cell_type_x', 'cell_id'], axis = 1)
    y_train = expression['cell_type_y']

# ## ML Pipeline start
RSEED = 42

steps = [('Scale', StandardScaler()), ('pca', PCA()), ('RandomForest', RandomForestClassifier())]
pipeline = Pipeline(steps)

if snakemake.wildcards['modality'] == 'scRNASeq' or snakemake.wildcards['modality'] == 'snRNASeq':
    parameters = {
        'pca__n_components': [25, 50],
        'RandomForest__max_features': [4, 6, 10],
        'RandomForest__n_estimators': [100, 150]
    }
else:
    parameters = {
        'pca__n_components': [20, 39],
        'RandomForest__max_features': [4, 6, 10],
        'RandomForest__n_estimators': [100, 150]
    }

print("starting grid search")
grid = GridSearchCV(pipeline, param_grid=parameters, cv=5, return_train_score=True)

grid.fit(X_train, y_train)

joblib.dump(grid.best_estimator_, snakemake.output['model'], compress = 1)

print(f"Training score: {grid.score(X_train, y_train)}")
print(f"The best parameters are: {grid.best_params_}")