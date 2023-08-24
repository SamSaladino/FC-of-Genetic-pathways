import numpy as np
import pandas as pd

""" This Module provide with useful methods to
transform gene expression data into pathway expression data
"""


def pathway_query(path,ensemble_ids):
    """ This function filter all pathways and that match with 
    the samples genes.

    Parameters:
    -----------
    path : str
        the excel file path with the corresponding ensemble IDs
        to the pathways
    ensemble_ids : file, list, dict, pandas DataFrame
        the samples ensemble IDs
    
    Returns:
    ---------
    pathway_query : The matching genes pathway DataFrame
    """
    # Read the pathway file from Reactome
    pathways = pd.read_excel(path, header=None)
    # Rename columns
    pathways.columns = ['Ensemble', 'Pathways']

    # Get the genes in corresponding pathways
    genes_annotated = pd.unique(pathways.Ensemble)
    # Get the matching genes
    genes_path = set(ensemble_ids) & set(genes_annotated)

    return pathways.query('Ensemble in @genes_path')

def pathway_expression(expression_matrix,pathway_names, pathway_data):
    """Transform the gene expression matrix into the
    pathway expression matrix

    Parameters:
    -----------
    expression_matrix : pandas DataFrame, numpy array
        gene expression matrix
    pathway_names : pandas DataFrame
        pathway names just that
    pathway_data : pandas DataFrame
        previous function executed

    Returns:
    --------
    pathway_expression : pandas DataFrame
        the pathway expression matrix
    """
    X = []
    # Run for each pathway
    for p in pathway_names:
        queried_pathways = pathway_data[pathway_data.Pathways == p]
        # Genes in this pathway
        g_pathways = queried_pathways.Ensemble
        # Calculate the pathway Expression
        matrix_values = np.sum(abs(expression_matrix[g_pathways]), axis=1) / g_pathways.shape[0]
        X.append(matrix_values)
    
    return pd.DataFrame(data=X, index=pathway_names)


if __name__ == '__main__':
    pass
