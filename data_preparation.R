setwd('C:/Users/zelda/mgr')


BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133plus2.db")
BiocManager::install("GSVA", update = TRUE)
BiocManager::install("GSEABase")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("KEGGREST")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GSEABenchmarkeR")
BiocManager::install("DataCombine")
BiocManager::install("genefilter")

remotes::install_version("RSQLite", version = "2.2.5")

library("DataCombine")
library("KEGGREST")
library("GSEABenchmarkeR")
library("EnrichmentBrowser")
library("GSVA")
library("GSEABase")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("hgu133a.db")    ##for Human
library("hgu133plus2.db")
library("plyr")
library("HDF5Array")
library("genefilter")
BiocManager::install("pd.hg.u133a")
library("pd.hg.u133a")
BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
options(connectionObserver = NULL)

## load data

geo2kegg = loadEData("geo2kegg")
## 42 datasety
## assayData: x features, y samples 
## - hat contains expression levels of x probe sets measured for y patients.


## summarize expression levels for probes annotated to the same gene

geo2kegg_sum = maPreproc(geo2kegg)
names(geo2kegg)


single_dataset = geo2kegg[[1]]




# Functional enrichment


filtered_eset <- nsFilter(single_dataset, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE,feature.exclude="^AFFX")
head(pData(single_dataset))
table(single_dataset$Group)

dataset_1_filtered = filtered_eset$eset

dataset_1_es = gsva(dataset_1_filtered, gsc,
                                  kcdf=c("Gaussian"),
                                  abs.ranking=FALSE,
                                  min.sz=2,
                                  max.sz=Inf, ## all paths even if short
                                  parallel.sz=1L,
                                  mx.diff=TRUE,
                                  ssgsea.norm=TRUE,
                                  verbose=FALSE)
