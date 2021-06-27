# Analysis code for Di Pilato et al., 2021 [1]

## Data
Single-cell transcriptome data used by the code in this repository are stored on GEO: [to be updated after GEO submission is processed by GEO]. Pan-cancer TCGA data from [2] can be found at https://gdc.cancer.gov/about-data/publications/pancanatlas. KP1.9 mouse lung adenocacinoma scRNA-seq data from [3] was obtained from GSE127465: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465 

## Methods
### Contributors 
* James Billingsley (JB)
* Marius Messemaker (MM)

### Single-cell transcriptome analysis 
#### From sequencing reads to raw gene counts (Contributor: JB)
Reads were processed using the inDrop v3 pipeline implemented in bcbio-nextgen version 1.2.4-76d5c4b (https://github.com/bcbio/bcbio-nextgen). Briefly, dual cell barcodes, sample barcodes and UMIs were detected using umis [4] (correcting cellular and sample barcodes of edit distance 1 from expected cell barcodes). Cell barcodes with fewer than 500 total reads assigned to them were discarded. Reads were assigned to Mus Musculus transcripts (GRCm38 (mm10) release M23) using Rapmap [5] and counts of reads per transcript per unique UMI were generated for each cell. 

#### Preprocessing (Contributor: JB)

#### Cell annotation (Contributor: MM)

#### Analyses using annotated data (Contributor: MM) 

#### Analysis of KP1.9 mouse lung adenocarcinoma scRNA-seq data [4] (Contributor: MM)

### Survival analysis (Contributor: MM) 

## References
[1] To be updated after publication.  
[2] Liu, J., Lichtenberg, T., Hoadley, K.A., Poisson, L.M., Lazar, A.J., Cherniack, A.D., Kovatich, A.J., Benz, C.C., Levine, D.A., Lee, A. V, et al. (2018). An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell 173, 400-416.e11.  
[3] Zilionis, R., Engblom, C., Pfirschke, C., Savova, V., Zemmour, D., Saatcioglu, H.D., Krishnan, I., Maroni, G., Meyerovitz, C. V, Kerwin, C.M., et al. (2019). Single-Cell Transcriptomics of Human and Mouse Lung Cancers Reveals Conserved Myeloid Populations across Individuals and Species. Immunity 50, 1317-1334.e10.  
[4] Svensson, V., Natarajan, K., Ly, LH. et al. (2017). Power analysis of single-cell RNA-sequencing experiments. Nat Methods 14, 381–387.  
[5] Avi Srivastava, Hirak Sarkar, Nitish Gupta, Rob Patro. (2016). RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes, Bioinformatics, Volume 32, Issue 12, Pages i192–i200.

