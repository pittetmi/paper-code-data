# Analysis code for Di Pilato et al., 2021 [1]

## Data
Single-cell transcriptome data used by the code in this repository are stored on GEO: [to be updated after GEO submission is processed by GEO]. Pan-cancer TCGA data from [2] can be found at https://gdc.cancer.gov/about-data/publications/pancanatlas. KP1.9 mouse lung adenocacinoma scRNA-seq data from [3] was obtained from GSE127465: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465 

## Methods
### Contributors 
* James Billingsley (JB)
* Marius Messemaker (MM)

### Single-cell transcriptome analysis 

| Mehtod description| Notebook | Contributor |
|     ---      |      ---       |     :---:     |
|Reads were processed using the inDrop v3 pipeline implemented in bcbio-nextgen version 1.2.4-76d5c4b (https://github.com/bcbio/bcbio-nextgen). Briefly, dual cell barcodes, sample barcodes and UMIs were detected using umis [4] (correcting cellular and sample barcodes of edit distance 1 from expected cell barcodes). Cell barcodes with fewer than 500 total reads assigned to them were discarded. Reads were assigned to Mus Musculus transcripts (GRCm38 (mm10) release M23) using Rapmap [5] and counts of reads per transcript per unique UMI were generated for each cell. |NA| JB |
|In this notebook the raw gene matrix is filtered and cell counts were normalized with the SCTransform function from the Seurat package [6]. To reduce dimensions, Principal Component Analysis was performed using the function RunPCA from the Seurat package set to 15 Principal Components. A nearest neighbor graph and UMAP were created with the functions FindNeighbors and runUMAP from the Seurat package.|[to be updated with the R markdown file]| JB |
|In this notebook we load the allntova (All cells, Not aPD1 Treated, pOVA immunogenic tumor) dataset from Seurat. In addition, we prepare reference scRNA-seq datasets that are used to annotate and validate cell state annotation of the allntova dataset. Reference datasets are linked in the notebook.|[part1b_prepare_annotation.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/part1b_prepare_annotation.ipynb)|MM|
|In this notebook cells in the allntova dataset are assigned to prior annotated reference states from multiple reference scRNA-seq datasets. In addition, single-cell expression of known marker genes from different cell states are plotted (Fig. 1C, 4A, 4C, 5A, 6E, S1D, & S5A). In this way, cell states can be identified in the UMAP from the allntova dataset and annotated by unsupervised clustering in part 2b annotation. Reference datasets are linked in the notebook. |[part2a_annotation.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/part2a_annotation.ipynb)|MM|
|In this notebook cells in the allntova dataset are clustered with unsupervised Leiden clustering at different resolutions. Next, appriopiate cluster resolutions are used for different cell populations (e.g. T cells) to acquire a final cell state annotation that matches classification by prior-annotated reference cell states and cell state marker genes (see part 2a).|[part2b_annotation.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/part2b_annotation.ipynb)|MM|
|In this notebook the final annotation of the T and Dendritic (DC) cell populations in the allntova dataset (see part 2b) is validated by comparing each final annotated cell state to each prior-annotated reference state in the datasets that were prepared in part 1b using a reciprocal similarity score [6] (Fig. S1B & S1C). In addition, we examine whether cell state enriched gene sets represent/match cell state annotation (Fig. S1E). Lastly, to make it possible for other researchers to review the single-cell transcriptome data, an interactive SPRING explorer is exported that will be loaded onto the SPRING server: https://github.com/AllonKleinLab/SPRING_dev|[part3_annotation_validation.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/part3_annotation_validation.ipynb)|MM|
|In this notebook we generate Fig. 1A, 1D, 4B & 6F using the annotated single-cell transcriptome data.|[part4_figures.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/part4_figures.ipynb)|MM|
|In this notebook we re-analyze KP1.9 mouse lung adenocarcinoma single-cell transcriptome data from [3] and obtained from GSE127465: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465. Fig. S4A-D, S4G, S5B, & S6B-C were generated using the code in this notebook.|[KP1_9.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/KP1_9.ipynb)|MM|

### Survival analysis 

| Method description | Notebook| Contributor |
| ------------- | ------------- | :---: |
| In this notebook we examine whether CXCR6 expression predicts survival in human cancer patients using bulk RNA-Seq data sets publically available from the TCGA pan-cancer dataset: https://gdc.cancer.gov/about-data/publications/pancanatlas. Fig. 7 & S7 were generated using the code in this notebook.   | [survival_analysis.ipynb](https://github.com/pittetmi/paper-code-data/blob/main/Di_Pilato_et_al_2021/survival_analysis.ipynb)  | MM |




## References
[1] To be updated after publication.  
[2] Liu, J., Lichtenberg, T., Hoadley, K.A., Poisson, L.M., Lazar, A.J., Cherniack, A.D., Kovatich, A.J., Benz, C.C., Levine, D.A., Lee, A. V, et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell.   
[3] Zilionis, R., Engblom, C., Pfirschke, C., Savova, V., Zemmour, D., Saatcioglu, H.D., Krishnan, I., Maroni, G., Meyerovitz, C. V, Kerwin, C.M., et al. (2019). Single-Cell Transcriptomics of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations across Individuals and Species. Immunity.   
[4] Svensson, V., Natarajan, K., Ly, LH. et al. (2017). Power analysis of single-cell RNA-sequencing experiments. Nat Methods.  
[5] Avi Srivastava, Hirak Sarkar, Nitish Gupta, Rob Patro. (2016). RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes, Bioinformatics.  
[6] Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome biology.  
[7] Gerhard, G.M., Bill, R., Messemaker, M., Klein, A.M., Pittet, M.J.n (2021). Tumor-infiltrating dendritic cell states are conserved across solid human cancers. J Exp Med.


