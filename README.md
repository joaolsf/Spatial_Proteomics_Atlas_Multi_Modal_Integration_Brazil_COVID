Brazil_COVID_cell_atlas
Longitudinal and full autopsy study to generate a spatially resolved single cell atlas of the post-mortem lung in COVID-19 and to integrate with ante-mortem signatures from clinical and high-plex peripheral blood profiling data. Here, you can find the codebase of the study.

# Disease trajectories in hospitalized COVID-19 patients are predicted by clinical and peripheral blood signatures representing distinct lung pathologies.
João Da Silva Filho1,2# , Vanessa Herder3# , Matthew P. Gibbins1,2# , Monique Freire dos Reis4.5# , Gisely Cardoso Melo6# , Michael J. Haley9# , Carla Cristina Judice10 , Fernando Fonseca Almeida Val5,6 , Mayla Borba5,13 , Tatyana Almeida Tavella10,15 , Vanderson de Sousa Sampaio6,14 , Charalampos Attipa1,11,16 , Fiona McMonagle1,12 , Marcus Vinicius Guimaraes de Lacerda6,7,8* , Fabio Trindade Maranhão Costa10* , Kevin N. Couper9* , Wuelton Marcelo Monteiro5,6* , Luiz Carlos de Lima Ferreira5,6* , Christopher Alan Moxon1*, Massimo Palmarini3*, Matthias Marti1,2*

# Check out our pre-print on medRxiv:
[doi: https://doi.org/10.1101/2023.09.08.23295024]

# Background
COVID-19 is characterized by a broad range of symptoms and disease trajectories. Understanding the correlation between clinical biomarkers and lung pathology over the course of acute COVID-19 is necessary to understand its diverse pathogenesis and inform more precise and effective treatments. Here, we present an integrated analysis of longitudinal clinical parameters, peripheral blood biomarkers, and lung pathology in COVID-19 patients from the Brazilian Amazon. We identified core clinical and peripheral blood signatures differentiating disease progression between recovered patients from severe disease and fatal cases. Signatures were heterogenous among fatal cases yet clustered into two patient groups: “early death” (< 15 days of disease until death) and “late death” (> 15 days). Progression to early death was characterized systemically and in lung histopathology by rapid, intense endothelial and myeloid activation/chemoattraction and presence of thrombi, driven by SARS-CoV-2 + macrophages. In contrast, progression to late death was associated with fibrosis, apoptosis and abundant SARS-CoV-2 + epithelial cells in post-mortem lung, with cytotoxicity, interferon and Th17 signatures only detectable in the peripheral blood 2 weeks into hospitalization. Progression to recovery was associated with systemic pro-lymphogenic, Th2 and anti-inflammatory-mediated responses. By integrating ante-mortem longitudinal systemic and spatial single-cell lung signatures, we defined an enhanced set of prognostic clinical parameters predicting disease outcome for guiding more precise and optimal treatments. To our knowledge, this is the first study of any acute respiratory infection, where serial clinical data and peripheral blood samples have been linked to histopathological and spatially-resolved single-cell investigations of post-mortem lung samples.

# Code Layout
## Analysis:
### 1- Scripts containing all analysis based on longitudinal ante-mortem clinical and high-plex peripheral blood data:
PCA, UMAP, hierarchical clustering, K-means clustering, EFA, trajectory inference, machine learning predictive models.
### 2- Scripts containing all imaging mass cytometry analysis of the post-mortem lungs:
Pre-processing, cell segmentation, extraction of single-cell features, single-cell and spatial statistics analysis.
### 3- Multi-modal data integration:
Multi-factor analysis (MFA), Multi-Omics Factor Analysis (MOFA), Method for the Functional Integration of Spatial and Temporal Omics data (MEFISTO) and Sparse Decomposition of Arrays (SDA).
## Figures:
Scripts containing all steps to regenerate figures from the study.



