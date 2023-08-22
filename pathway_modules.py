import numpy as np
import pandas as pd
from scipy.stats import gmean


# class PathwayMatrix:
#     def __init__(self, folder_path,ensemble_id):
#         self.folder_path = folder_path
#         self.ensemble_id = ensemble_id


#     def transform(self,df):
#         self.matrix = []

#         for p in pathways_names:


#     def _pathway_df(self):
#         pathways = pd.read_excel(self.folder_path,header=None)
#         pathways.columns = ["Ensemble","Pathways"]


pathways_names = pd.unique(pathways.Pathways)

genes_annotated = pd.unique(pathways.Ensemble)

# Genes in pathways
genes_path = set(ensemble_id) & set(genes_annotated)

# Pathways with annotated genes
pathways_annotated = top_pathways.query('Ensemble in @genes_path')

def Matrix(col=np.arange()):

X = []
for pathway in pathways_names:
    # genes evaluated in pathway (variable to operate)
    g_pathways = df[
        pathways_annotated[pathways_annotated.Pathways == pathway].Ensemble
    ]
    # Ecuation above
    X.append(np.sum(abs(g_pathways),axis=1)/g_pathways.shape[1])

X = pd.DataFrame(data=X, columns=tissues,index=pathways_names)