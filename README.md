# Testing-on-BayesPrism
Project for STAT 556

## File information
excute_pancreas.r: R script for generate cell type fraction on BayesPrism package for pancreas dataset.

excute_PBMC.r: R script for generate cell type fraction on BayesPrism package for PBMC dataset.

bulk_RNA-seq_deconvolution_BayesPrism_CybersortX.Rmd: R markdown for data analysis, detail information included.

bulk_RNA-seq_deconvolution_BayesPrism_CybersortX.html: interpreted R markdown in html format (for data analysis, detail information included). 

## Important dataset:

#### Notation:
B: bulk gene expression matrix, represents expression levels of the g-th gene of the sth sample.

A: signifture matrix, g-th reference gene expression level for the k-th cell type.

P: unobservable truth of cell type proportions in the samples

ref: reference count matrix file with gene name(row) * cell id (col)


#### 1.PBMC
B.PBMC:    "PBMC_ARCHIVE/Fig2b-WholeBlood_RNAseq.txt"

A.PBMC:    "PBMC_ARCHIVE/Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt"

ref.PBMC:  "PBMC_ARCHIVE/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt"

cyb_P.PBMC:"PBMC_ARCHIVE/CIBERSORTx_PBMC_Adjusted.txt" -- obtained in cybersortX website

bp_P.PBMC: "PBMC_ARCHIVE/PBMC_theta.rds" -- obtained in r package (BayesPrism)

true_P.PBMC:"PBMC_ARCHIVE/ground_truth_whole_blood.txt"

#### 2.pancreas
B.pancreas:"pancreas_ARCHIVE/pancreas_bulk_v2.txt"

A.pancreas: "pancreas_ARCHIVE/pancreas_inferred_phenoclasses.txt" -- obtained in cybersortX website

ref.pancreas:"pancreas_ARCHIVE/pancreas_scref_v2.txt"

cyb_P.pancreas:"pancreas_ARCHIVE/CIBERSORTx_pancreas_Adjusted.txt" -- obtained in cybersortX website

bp_P.pancreas:"pancreas_ARCHIVE/pancreas_theta.rds" -- obtained in r package (BayesPrism)

true_P.pancreas:"pancreas_ARCHIVE/pancreas_truth.txt"

## Procedure for replicate the result:
step 1. run pancreas_excute.r to obtain "bp_P.pancreas"

step 2. run PBMC_excute.r to obtain "bp_P.PBMC"

step 3. run "bulk_RNA-seq_deconvolution_BayesPrism_CybersortX"

Note that step 1 & 2 can be skipped, I have included the data in the file already.
