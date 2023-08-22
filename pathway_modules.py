import numpy as np
import pandas as pd
from scipy.stats import gmean


class PathwayMatrix:
    def __init__(self) -> None:
        

# Read the Top 28 pathways excel list
top_pathways = pd.read_excel("./data_external/data_pathways/pathways2_sorted.xls",header=None)
# Name the columns
top_pathways.columns = ["Ensemble","Pathways"]

# Count pathways
pathways_names = pd.unique(top_pathways.Pathways)
# Count genes in pathways
genes_annotated = pd.unique(top_pathways.Ensemble)