""" This Module provides useful methods to
transform gene expression data into pathway expression data

For more details on this see author(s):
Costa S., Nieves J.
"""

import numpy as np
import pandas as pd

from typing import Optional, List
from scipy.stats import gmean


class PathwayAnalysis:
    """
    A class to perform pathway analysis using pathway and gene expression data.
    Attributes
    ----------
    pathway_file : str
        Path to the file containing pathway data.
    expression_file : str
        Path to the file containing gene expression data.
    Methods
    -------
    pathway_query(gene_ids: list) -> pd.DataFrame
        Queries the pathways DataFrame for the given gene IDs and returns
        the matching pathways.
    pathway_expression(pathway_names: list, pathway_data: pd.DataFrame) -> pd.DataFrame
        Calculates the expression values for the given pathway names using the
        provided pathway data.
    """
    
    def __init__(self, pathways: Optional[pd.DataFrame] = None,
                 pathway_file: Optional[str] = None):
        """
        Initializes the PathwayAnalysis class with the given pathway and
        expression files.
        Parameters
        ----------
        pathway_file : str
            Path to the file containing pathway data.
        """
        if pathways is not None:
            self._pathways: pd.DataFrame = pathways.copy()
        elif pathway_file is not None:
            self._pathways: pd.DataFrame = pd.read_csv(pathway_file, names=['Gene', 'Pathway'])
        else:
            raise TypeError('Pathways data or pathway file is needed')

        self.pathways: pd.DataFrame = self._pathways.copy()
        self.genes_annotated: np.ndarray = self.pathways.Gene.unique()

    def pathway_query(self, gene_ids: List[str]) -> pd.DataFrame:
        """
        Queries the pathways DataFrame for the given gene IDs and returns
        the matching pathways.
        Parameters
        ----------
        gene_ids : list
            List of gene IDs to query.
        Returns
        -------
        pd.DataFrame
            DataFrame containing the pathways that match the given gene IDs.
        """
        common_genes: List[str] = list(set(gene_ids) & set(self.genes_annotated))
        self.pathways = self._pathways.query('Gene in @common_genes').copy()
        self.common_genes: np.ndarray = np.array(common_genes)
        self.pathway_names = list(set(self.pathway_names) & set(self.pathways.Pathway.unique()))
        self.pathway_names.sort()
        return self.pathways

    def calc_fold_change(self, 
                         gene_expression: Optional[pd.DataFrame] = None,
                         gene_expression_file: Optional[str] = None,
                         metadata: Optional[pd.DataFrame] = None,
                         metadata_file: Optional[str] = None,
                         normal_label: str = '0',) -> None:
        '''
        gene_expression: Optional[pd.DataFrame]
            Dataframe containing the gene expression data
        expression_file: Optional[str]
            Path to the file containing the gene expression data.
        metadata: Optional[pd.DataFrame]
            Dataframe containing the the label which each sample belong
        metadata_file: Optional[str]
            Path to the file containing the label which each sample belong
        normal_label: Optional[str]
            Label of normal group
        '''
        if gene_expression is not None:
            self.gene_expression = gene_expression.copy()
        elif gene_expression_file is not None:
            self.gene_expression = pd.read_csv(gene_expression_file, index_col=0)
        else:
            raise TypeError('The dataframe or the file path with the gene expression data is needed')

        if metadata is not None:
            self.metadata = metadata.copy()
        elif metadata_file is not None:
            self.metadata = pd.read_csv(metadata_file, dtype=str)
        else:
            raise TypeError('The dataframe or the file path with the metadata is needed')


        normal_mask = self.metadata.iloc[:, 1] == normal_label
        normal_sample = self.metadata[normal_mask].iloc[:, 0]
        normal_data = self.gene_expression[normal_sample] + 0.1

        disease_mask = ~normal_mask
        disease_sample = self.metadata[disease_mask].iloc[:, 0]
        disease_data = self.gene_expression[disease_sample] + 0.1

        ref = gmean(normal_data, axis=1)
        self.fold_change = np.log2(disease_data.div(ref, axis=0))

    def pathway_expression(self,
                           gene_expression: Optional[pd.DataFrame] = None,
                           gene_expression_file: Optional[str] = None,
                           metadata: Optional[pd.DataFrame] = None,
                           metadata_file: Optional[str] = None,
                           normal_label: str = '0',
                           pathway_names: Optional[List[str]] = None
                           ) -> pd.DataFrame:
        """
        Calculates the expression values for the given pathway names using the
        provided pathway data.
        Parameters
        ----------
        gene_expression: Optional[pd.DataFrame]
            Dataframe containing the gene expression data
        expression_file : Optional[str]
            Path to the file containing the gene expression data.
        metadata: Optional[pd.DataFrame]
            Dataframe containing the the label which each sample belong
        metadata_file: Optional[str]
            Path to the file containing the label which each sample belong
        normal_label: Optional[str]
            Label of normal group
        pathway_names : Optional[list]
            List of pathway names to calculate expression values for. If it is 
            not provided all pathways in are used
        Returns
        -------
        pd.DataFrame
            DataFrame containing the calculated expression values for each
            pathway.
        """
        if pathway_names is None:
            pathway_names = self.pathways.Pathway.unique()
            pathway_names.sort()
        self.pathway_names: List[str] = pathway_names

        if not hasattr(self, 'fold_change') or self.fold_change is None:
            self.calc_fold_change(
                gene_expression=gene_expression,
                gene_expression_file=gene_expression_file,
                metadata=metadata,
                metadata_file=metadata_file,
                normal_label=normal_label,
            )

        gene_ids: np.ndarray = self.fold_change.index.unique().to_numpy()
        self.pathway_query(gene_ids)
        X: np.ndarray = np.zeros((len(self.pathway_names), self.fold_change.shape[1]))
        for i, pname in enumerate(self.pathway_names):
            queried_pathway: pd.DataFrame = self.pathways[self.pathways.Pathway == pname]
            g_pathway: pd.Series = queried_pathway.Gene
            X[i] = self.fold_change.loc[g_pathway].abs().mean()

        self.pathway_data: pd.DataFrame = pd.DataFrame(
            data=X, index=self.pathway_names, columns=self.fold_change.columns)

        return self.pathway_data


if __name__ == '__main__':
    pass
