library(snowfall)
library(NMF)
library(gplots)
library(scran)
library(BiocParallel)
library(devtools)
suppressWarnings(library(BayesPrism))
PBMCbulk <- read.table(file = "PBMC_ARCHIVE/Fig2b-WholeBlood_RNAseq.txt", header = TRUE,sep = "\t")
rownames(PBMCbulk) <- PBMCbulk[,1]
PBMCbulk<-PBMCbulk[,-1]

PBMCref <- read.table(file = "PBMC_ARCHIVE/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt", header = TRUE,sep = "\t")
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
bp.res_PBMC
# extract posterior mean of cell type fraction theta
theta_PBMC <- get.fraction (bp=bp.res_PBMC,
                            which.theta="final",
                            state.or.type="type")

print(head(theta_PBMC))

saveRDS(theta_PBMC, "PBMC_ARCHIVE/PBMC_theta.rds")

