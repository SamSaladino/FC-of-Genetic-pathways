import numpy as np
import pandas as pd



class PathwayAnalyzer:

    def __init__(self,path):
        self.pathways = pd.read_excel(path,header=None)
        self.pathways.columns = ['Pathways','Ensemble']

# # class PathwayMatrix:
# #     def __init__(self, folder_path,ensemble_id):
# #         self.folder_path = folder_path
# #         self.ensemble_id = ensemble_id


# #     def transform(self,df):
# #         self.matrix = []

# #         for p in pathways_names:


# #     def _pathway_df(self):
# #         pathways = pd.read_excel(self.folder_path,header=None)
# #         pathways.columns = ["Ensemble","Pathways"]


# def pathways_query(path=str, ensemble_id=list):
#     pathways = pd.read_excel(path, header=None)
#     pathways.columns = ["Ensemble","Pathways"]
#     genes_annotated = pd.unique(pathways.Ensemble)
#     genes_path = set(ensemble_id) & set(genes_annotated)
    
#     return pathways.query('Ensemble in @genes_path')



# def matrix(df,col=False, col_name=None):
#     pathways_names = pd.unique(pathways_query.Pathways)
#     X = []
#     for p in pathways_names:
#         g_pathways = df[
#             pathways_query[
#                 pathways_query.Pathways == p
#                 ].Ensemble
#             ]
#         X.append(np.sum(abs(g_pathways),axis=1)/g_pathways.shape[1])
    
#     if col:
#         return pd.DataFrame(data=X, columns=col_name,index=pathways_names)
#     else:
#         return pd.DataFrame(data=X,index=pathways_names)
    import numpy as np
import pandas as pd

class PathwayAnalyzer:
    def __init__(self, path):
        self.pathways = pd.read_excel(path, header=None)
        self.pathways.columns = ["Ensemble", "Pathways"]
        self.pathway_names = pd.unique(self.pathways.Pathways)
        
    def pathways_query(self, ensemble_ids):
        genes_annotated = pd.unique(self.pathways.Ensemble)
        genes_path = set(ensemble_ids) & set(genes_annotated)
        return self.pathways.query('Ensemble in @genes_path')
    
    def calculate_matrix(self, data_frame, col=False, col_name=None):
        X = []
        for p in self.pathway_names:
            g_pathways = self.pathways_query(self.pathways[self.pathways.Pathways == p].Ensemble)
            # Modify the following calculation as needed
            matrix_values = np.sum(abs(data_frame[g_pathways]), axis=1) / g_pathways.shape[1]
            X.append(matrix_values)
        
        if col:
            return pd.DataFrame(data=X, columns=col_name, index=self.pathway_names)
        else:
            return pd.DataFrame(data=X, index=self.pathway_names)

# Usage
pathway_analyzer = PathwayAnalyzer("path_to_excel_file.xlsx")
ensemble_ids = [list_of_ensemble_ids]
data_frame = pd.DataFrame(...)  # Provide your data frame here
result_matrix = pathway_analyzer.calculate_matrix(data_frame, col=True, col_name=["Pathway1", "Pathway2", ...])
print(result_matrix)
