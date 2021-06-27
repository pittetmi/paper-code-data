# Analysis code for Di Pilato et al., 2021 [1]

## Data
Single-cell transcriptome data used by the code in this repository are stored on GEO. Publically available pan-cancer TCGA data from [2] can be found at https://gdc.cancer.gov/about-data/publications/pancanatlas. Publically available KP1.9 mouse lung adenocacinoma scRNA-seq data was obtained from GEO

## Methods
### Contributors 
* James Billingsley (JB)
* Marius Messemaker (MM)

### Single-cell transcriptome analysis 
#### From sequencing reads to raw gene counts (Contributor: JB)
Reads were processed using the inDrop v3 pipeline implemented in bcbio-nextgen version 1.2.4-76d5c4b (https://github.com/bcbio/bcbio-nextgen). Briefly, dual cell barcodes, sample barcodes and UMIs were detected using umis [3] (correcting cellular and sample barcodes of edit distance 1 from expected cell barcodes). Cell barcodes with fewer than 500 total reads assigned to them were discarded. Reads were assigned to Mus Musculus transcripts (GRCm38 (mm10) release M23) using Rapmap [3] and counts of reads per transcript per unique UMI were generated for each cell. 

#### Preprocessing (Contributor: JB)

#### Cell annotation (Contributor: MM)

#### Analyses using annotated data (Contributor: MM) 

#### Analysis of KP1.9 mouse lung adenocarcinoma scRNA-seq data [4] (Contributor: MM)

### Survival analysis (Contributor: MM) 

## References
