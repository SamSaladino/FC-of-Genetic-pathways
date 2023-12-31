# FC-of-Genetic-pathways
### Log Fold Change of Genetic pathways, Top Pathways in cancer TCGA:

This is the official repository associated to the paper [**Average fold-change of genetic pathways
in biological transitions**](https://doi.org/10.48550/arXiv.2305.11245)

We present a quantitative metric aimed at assessing the average fold-change within genetic pathways during biological transitions. This metric was employed to broadly characterize the shift from a healthy tissue to a tumor state.
We use gene expression data from the TCGA portal for 16 tumors and the
Reactome compilation of pathways.

### Notebooks and scripts:
 In this folder you can find two python notebooks and one python script

 * *<span style="color:red">Expr_pathways.ipynb</span>* Generate the figures and tables corresponding to all 16 tumor samples using the top 28 pathways from Reactome database.

 * *<span style="color:red">Melanoma_exp_pathways.ipynb</span>* Generate the data means, tables and graphics corresponding to melanoma data.
  
 * *<span style="color:red">pathway_modules.py</span>* is a module imported in the notebooks for easier manipulation

### Folders
* *<span style="color:blue">data_external</span>* : contains the pathways reactome info (top 28 pathways and all pathways) and the melanoma data downloaded from TCGA

* *<span style="color:blue">data_generated</span>* : data generated in this repository or DarioALeonValido/evolp repository
  
* *<span style="color:blue">figures_tables</span>* : figures and tables of the paper generated in the notebooks


*I hope this formatting works for you! If you have any more questions or need further assistance, feel free to ask.*