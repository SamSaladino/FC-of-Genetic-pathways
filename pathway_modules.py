import numpy as np
import pandas as pd


#################################################
def pathway_read(path):
    return pd.read_excel(path, header=None)

def genes_coincidences(path, ensemble_ids):
    pathways = pd.read_excel(path, header=None)
    pathways.columns = ['Pathways', 'Ensemble']

    genes_annotated = pd.unique(pathways.Ensemble)
    return set(ensemble_ids) & set(genes_annotated)



########################################################
def pathway_query(path,ensemble_ids):
    pathways = pd.read_excel(path, header=None)
    pathways.columns = ['Pathways', 'Ensemble']

    genes_annotated = pd.unique(pathways.Ensemble)
    genes_path = set(ensemble_ids) & set(genes_annotated)

    return pathways.query('Ensemble in @genes_path')

def pathway_expression(expression_matrix,pathway_names, pathway_data):
    X = []
    for p in pathway_names:
        queried_pathways = pathway_data[pathway_data.Pathways == p]
        g_pathways = queried_pathways.Ensemble
        matrix_values = np.sum(abs(expression_matrix[g_pathways]), axis=1) / g_pathways.shape[0]
        X.append(matrix_values)
    
    return pd.DataFrame(data=X, index=pathway_names)


if __name__ == '__main__':
    print("run")
