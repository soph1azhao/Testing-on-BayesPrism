install.packages("NMF")
install.packages("gplots")
install.packages("scran")
install.packages("BiocParallel")
install.packages("devtools")
library(snowfall)
library(NMF)
library(gplots)
library(scran)
library(BiocParallel)
library(devtools)
install_github("Danko-Lab/BayesPrism/BayesPrism")
suppressWarnings(library(BayesPrism))

setwd("S:/Study/uvic/Fall 2023/STAT-556-Directed-Study/BayesPrism-main")

##PBMC Preprocess##################################################################################
PBMCbulk <- read.table(file = "PBMC/Fig2b-WholeBlood_RNAseq.txt", header = TRUE,sep = "\t")
rownames(PBMCbulk) <- PBMCbulk[,1]
PBMCbulk<-PBMCbulk[,-1]

PBMCref <- read.table(file = "PBMC/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt", header = TRUE,sep = "\t")
rownames(PBMCref) <- PBMCref[,1]
PBMCref<-PBMCref[,-1]
head(PBMCref[,1:3])
head(PBMCbulk[,1:3])

bk.dat <- t(as.matrix(PBMCbulk,byrow = TRUE,
                       dimnames = list(rownames(PBMCbulk),
                                       colnames(PBMCbulk))))
#head(bk.dat[,1:5])
head(bk.dat[,1:5])
head(PBMCref[,1:5])
sc.dat <- t(as.matrix(PBMCref,byrow = TRUE,
                       dimnames = list(rownames(PBMCref),
                                       colnames(PBMCref))))
head(PBMCref[,1:5])


PBMCref_cell <- colnames(PBMCref)
# Cell types
PBMCcell_types <- c("B.cells", "T.cells.CD4", "T.cells.CD8", "NKT.cells", "NK.cells", "Monocytes")

library(stringr)
# Extract cell type and additional information using stringr
matches <- str_match(PBMCref_cell, paste0("(", paste(PBMCcell_types, collapse = "|"), ")\\.?(.*)"))
cell.type.labels<-matches[, 2]
cell.state.labels<-gsub("\\.", "-", matches[, 1])

sort(table(cell.type.labels))
sort(table(cell.state.labels))
table(cbind.data.frame(cell.state.labels, cell.type.labels))


##pancreas Preprocess##################################################################################
pancreasbulk <- read.table(file = "pancreas/pancreas_bulk_v2.txt", header = TRUE,sep = "\t")
rownames(pancreasbulk) <- pancreasbulk[,1]
pancreasbulk<-pancreasbulk[,-1]

pancreasref <- read.table(file = "pancreas/pancreas_scref_v2.txt", header = TRUE,sep = "\t")
rownames(pancreasref) <- pancreasref[,1]
pancreasref<-pancreasref[,-1]
head(pancreasref[,1:3])
head(pancreasbulk[,1:3])
bk.dat2 <- t(as.matrix(pancreasbulk,byrow = TRUE,
                       dimnames = list(rownames(pancreasbulk),
                                       colnames(pancreasbulk))))
#head(bk.dat[,1:5])
head(bk.dat2[,1:5])
head(pancreasref[,1:5])
sc.dat2 <- t(as.matrix(pancreasref,byrow = TRUE,
                       dimnames = list(rownames(pancreasref),
                                       colnames(pancreasref))))
#head(sc.dat[,1:5])
head(sc.dat2[,1:5])


pancreasref_cell <- colnames(pancreasref)
# Cell types
cell_types2 <- c("delta", "alpha", "epsilon", "ductal", "acinar", "beta", "PSC", "endothelial", "gamma","mast")

library(stringr)
# Extract cell type and additional information using stringr
matches2 <- str_match(pancreasref_cell, paste0("(", paste(cell_types2, collapse = "|"), ")\\.?(.*)"))
cell.type.labels2<-matches2[, 2]
cell.state.labels2<-gsub("\\.", "-", matches2[, 1])

# sort(table(cell.type.labels))
# sort(table(cell.state.labels))
# table(cbind.data.frame(cell.state.labels, cell.type.labels))

sort(table(cell.type.labels2))
sort(table(cell.state.labels2))
table(cbind.data.frame(cell.state.labels2, cell.type.labels2))

##PBMC outliers##################################################################################
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting.
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)
dim(sc.dat.filtered)

# #note this function only works for human data. For other species, you are advised to make plots by yourself.
# plot.bulk.vs.sc (sc.input = sc.dat.filtered2,
#                  bulk.input = bk.dat2
#                  #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
# )


sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")
##pancreas outliers##################################################################################
sc.stat2 <- plot.scRNA.outlier(
  input=sc.dat2, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels2,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting.
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)
bk.stat2 <- plot.bulk.outlier(
  bulk.input=bk.dat2,#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=sc.dat2, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels2,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
sc.dat.filtered2 <- cleanup.genes (input=sc.dat2,
                                   input.type="count.matrix",
                                   species="hs",
                                   gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                   exp.cells=5)
dim(sc.dat.filtered2)

# #note this function only works for human data. For other species, you are advised to make plots by yourself.
# plot.bulk.vs.sc (sc.input = sc.dat.filtered2,
#                  bulk.input = bk.dat2
#                  #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
# )


sc.dat.filtered.pc2 <-  select.gene.type (sc.dat.filtered2,
                                          gene.type = "protein_coding")
# # Select marker genes (Optional)
# # performing pair-wise t test for cell states from different cell types
#
# diff.exp.stat2 <- get.exp.stat(sc.dat=sc.dat2[,colSums(sc.dat2>0)>3],# filter genes to reduce memory use
#                                           cell.type.labels=cell.type.labels2,
#                                           cell.state.labels=cell.state.labels2,
#                                           pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
#                                           cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
#                                           n.cores=1 #number of threads
#                                           )
#
# sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc2,
#                                          stat=diff.exp.stat,
#                                          pval.max=0.01,
#                                          lfc.min=0.1)
# dim(sc.dat.filtered.pc.sig)


##PBMC runBPrism########################################################################################
myPrism_PBMC <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=bk.dat,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
bp.res_PBMC <- run.prism(prism = myPrism_PBMC, n.cores=50)
slotNames(bp.res_PBMC)


# extract posterior mean of cell type fraction theta
theta_PBMC <- get.fraction (bp=bp.res_PBMC,
                        which.theta="final",
                        state.or.type="type")

head(theta_PBMC)
# extract coefficient of variation (CV) of cell type fraction
theta.cv_PBMC <- bp.res_PBMC@posterior.theta_f@theta.cv

head(theta.cv_PBMC)

##pancreas runBPrism#######################################################################################
## Construct a prism object
myPrism2 <- new.prism(
  reference=sc.dat.filtered.pc2,
  mixture=bk.dat2,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels2,
  cell.state.labels = cell.state.labels2,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
bp.res2 <- run.prism(prism = myPrism2, n.cores=50)

slotNames(bp.res2)

# extract posterior mean of cell type fraction theta
theta2 <- get.fraction (bp=bp.res2,
                       which.theta="final",
                       state.or.type="type")

head(theta2)
# extract coefficient of variation (CV) of cell type fraction
theta.cv2 <- bp.res2@posterior.theta_f@theta.cv

head(theta.cv2)
# # extract posterior mean of cell type-specific gene expression count matrix Z
# Z.tumor <- get.exp (bp=bp.res2,
#                     state.or.type="type",
#                     cell.name="tumor")
#
# head(t(Z.tumor[1:5,]))
