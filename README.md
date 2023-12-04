# SurvIAE
In Diffuse Large B-Cell Lymphoma (DLBCL), several methodologies are emerging to derive novel biomarkers to
be incorporated in the risk assessment. We realized a pipeline that relies on autoencoders (AE) and 
Explainable Artificial Intelligence (XAI), to stratify prognosis and derive a gene-based signature.

AE was exploited to learn an unsupervised representation of the gene expression (GE) from three publicly available 
datasets, each with its own technology. Multi-layer perceptron (MLP) was used to classify prognosis from latent
representation. GE data were preprocessed as normalized, scaled, and standardized.
Four different AE architectures (Large, Medium, Small and Extra Small) were compared to find the most suitable
for GE data. The joint AE-MLP classified patients on six different outcomes: overall survival at 12, 36, 60 months
and progression-free survival (PFS) at 12, 36, 60 months. XAI techniques were used to derive a gene-based signature
aimed at refining the Revised International Prognostic Index (R-IPI) risk, which was validated in a fourth independent
publicly available dataset. 

We named our tool *SurvIAE: Survival prediction with Interpretable AE*.

If you find this repository useful for your research, please cite our paper in *Computer Methods and Programs in Biomedicine*:
```
@article{ZACCARIA2023107966,
title = {SurvIAE: Survival prediction with Interpretable Autoencoders from Diffuse Large B-Cells Lymphoma gene expression data},
journal = {Computer Methods and Programs in Biomedicine},
pages = {107966},
year = {2023},
issn = {0169-2607},
doi = {https://doi.org/10.1016/j.cmpb.2023.107966},
url = {https://www.sciencedirect.com/science/article/pii/S0169260723006326},
author = {Gian Maria Zaccaria and Nicola Altini and Giuseppe Mezzolla and Maria Carmela Vegliante and Marianna Stranieri and Susanna Anita Pappagallo and Sabino Ciavarella and Attilio Guarini and Vitoantonio Bevilacqua},
keywords = {Gene expression data, Survival prediction, Autoencoder, Explainable Artificial Intelligence}
}
```