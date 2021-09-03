

BiocManager::install(version = '3.12')

BiocManager::install("GSVA", version = "devel")
BiocManager::install("ggplot2")
require(devtools)
install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org",force = TRUE)
BiocManager::install("ggplot2", force = TRUE)
install.packages("remotes")
library(remotes)
install_github("salmonellix/GSVA", force=TRUE)
library("GSVA")
library("BiocParallel")
library("BiocSingular")

devtools::install_github(c("hadley/ggplot2", "GuangchuangYu/ggtree"), force = TRUE)


remotes::install_github("davismcc/scater",force = TRUE)
install.packages("devtools")
devtools::install_github("tidyverse/ggplot2")
library("ggplot2")
library("scater")
library("davismcc/scater")
install.packages("tidyverse")
BiocManager::install("scater", force = TRUE)
install.packages("umap")
library("umap")

install.packages("tsne")
library("tsne")
library("dplyr")
library("tidyr")
library("Seurat")

methods_gsva = c("plage_pca")








for (m in methods_gsva){
  
  for (i in 28:length(geo2kegg)){
    # read from file expression set: rows - enterezID, cols - genes
    file_name <- paste("Results/ds_expr_enterez2/ds_uniq", "_",i,"_e.csv", sep="")
    ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
    # main function
    ES_GSVA = gsva(data.matrix(ds_tmp), gsc,
                   method= c("plage_pca"),
                   kcdf=c("Gaussian"),
                   abs.ranking=FALSE,
                   min.sz=2,
                   max.sz=Inf, ## all paths even if short
                   parallel.sz=1L,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   verbose=FALSE)
    # write to apropriate folder
    file_out <- paste("Results/modifications/plage//ds", "_",i,"_zf.csv", sep="")
    write.csv(ES_GSVA, file_out, row.names = TRUE)
  }
}

#("gsva", "ssgsea", "zscore","zscore_stouffer","zscore_fisher", "plage", "plage_pca"),


ES_GSVA_pca = gsva(data.matrix(ds_tmp), gsc,
               method= c("plage_pca"),
               kcdf=c("Gaussian"),
               abs.ranking=FALSE,
               min.sz=2,
               max.sz=Inf, ## all paths even if short
               parallel.sz=1L,
               mx.diff=TRUE,
               ssgsea.norm=TRUE,
               verbose=FALSE)

ES_GSVA_umap = gsva(data.matrix(ds_tmp), gsc,
                   method= c("plage_umap"),
                   kcdf=c("Gaussian"),
                   abs.ranking=FALSE,
                   min.sz=2,
                   max.sz=Inf, ## all paths even if short
                   parallel.sz=1L,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   verbose=FALSE)

ES_GSVA_tsne = gsva(data.matrix(ds_tmp), gsc,
                   method= c("plage_tsne"),
                   kcdf=c("Gaussian"),
                   abs.ranking=FALSE,
                   min.sz=2,
                   max.sz=Inf, ## all paths even if short
                   parallel.sz=1L,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   verbose=FALSE)

ES_GSVA_zscore = gsva(data.matrix(ds_tmp), gsc,
                    method= c("zscore"),
                    kcdf=c("Gaussian"),
                    abs.ranking=FALSE,
                    min.sz=2,
                    max.sz=Inf, ## all paths even if short
                    parallel.sz=1L,
                    mx.diff=TRUE,
                    ssgsea.norm=TRUE,
                    verbose=FALSE)

ES_GSVA_gsva = gsva(data.matrix(ds_tmp), gsc,
                      method= c("gsva"),
                      kcdf=c("Gaussian"),
                      abs.ranking=FALSE,
                      min.sz=2,
                      max.sz=Inf, ## all paths even if short
                      parallel.sz=1L,
                      mx.diff=TRUE,
                      ssgsea.norm=TRUE,
                      verbose=FALSE)

ES_GSVA_ssgsea = gsva(data.matrix(ds_tmp), gsc,
                    method= c("ssgsea"),
                    kcdf=c("Gaussian"),
                    abs.ranking=FALSE,
                    min.sz=2,
                    max.sz=Inf, ## all paths even if short
                    parallel.sz=1L,
                    mx.diff=TRUE,
                    ssgsea.norm=TRUE,
                    verbose=FALSE)

ES_GSVA_stouffer = gsva(data.matrix(ds_tmp), gsc,
                      method= c("zscore_stouffer"),
                      kcdf=c("Gaussian"),
                      abs.ranking=FALSE,
                      min.sz=2,
                      max.sz=Inf, ## all paths even if short
                      parallel.sz=1L,
                      mx.diff=TRUE,
                      ssgsea.norm=TRUE,
                      verbose=FALSE)


ES_GSVA_fisher = gsva(data.matrix(ds_tmp), gsc,
                        method= c("zscore_fisher"),
                        kcdf=c("Gaussian"),
                        abs.ranking=FALSE,
                        min.sz=2,
                        max.sz=Inf, ## all paths even if short
                        parallel.sz=1L,
                        mx.diff=TRUE,
                        ssgsea.norm=TRUE,
                        verbose=FALSE)

ES_GSVA_plage = gsva(data.matrix(ds_tmp), gsc,
                      method= c("plage"),
                      kcdf=c("Gaussian"),
                      abs.ranking=FALSE,
                      min.sz=2,
                      max.sz=Inf, ## all paths even if short
                      parallel.sz=1L,
                      mx.diff=TRUE,
                      ssgsea.norm=TRUE,
                      verbose=FALSE)

pcavectorgset <- function(gSetIdx, Z) {
  if(is(Z, "dgCMatrix")){
    s <- BiocSingular::runPCA(Z[gSetIdx, ], rank=30)
  } else {
    s <- BiocSingular::runPCA(Z[gSetIdx, ], rank=30)
  }
  # first pca component
  s$rotation[ , 1]
}

rightsingularsvdvectorgset <- function(gSetIdx, Z) {
  if(is(Z, "dgCMatrix")){
    s <- BiocSingular::runExactSVD(Z[gSetIdx, ], rank=2)
  } else {
    s <- svd(Z[gSetIdx, ])
  }
  # first svd component
  s$v[, 1]
}



plage <- function(X, geneSets, parallel.sz, verbose=TRUE,
                  BPPARAM=SerialParam(progressbar=verbose)) {
  if(is(X, "dgCMatrix")){
    message("Please bear in mind that this method first scales the values of the gene
    expression data. In order to take advantage of the sparse Matrix type, the scaling
    will only be applied to the non-zero values of the data. This is a provisional 
    solution in order to give support to the dgCMatrix format.")
    
    Z <- Matrix::t(X)
    Z <- .dgCapply(Z, scale, 2)
    Z <- Matrix::t(Z)
    
    es <- bplapply(geneSets, pcavectorgset, Z,
                   BPPARAM=BPPARAM)
    
    es <- do.call(rbind, es)
    
    es <- as(es, "dgCMatrix")
    
  } else {
    
    Z <- t(scale(t(X)))
    
    es <- bplapply(geneSets, pcavectorgset, Z,
                   BPPARAM=BPPARAM)
    
    es <- do.call(rbind, es)
    
    if (length(geneSets) == 1)
      es <- matrix(es, nrow=1)
    
    rownames(es) <- names(geneSets)
    colnames(es) <- colnames(X)
  }
  
  es
}


plage_pca <- function(X, geneSets, parallel.sz, verbose=TRUE,
                      BPPARAM=SerialParam(progressbar=verbose)) {
  if(is(X, "dgCMatrix")){
    message("Please bear in mind that this method first scales the values of the gene
    expression data. In order to take advantage of the sparse Matrix type, the scaling
    will only be applied to the non-zero values of the data. This is a provisional 
    solution in order to give support to the dgCMatrix format.")
    
    Z <- Matrix::t(X)
    Z <- .dgCapply(Z, scale, 2)
    Z <- Matrix::t(Z)
    
    es <- bplapply(geneSets, pcavectorgset, Z,
                   BPPARAM=BPPARAM)
    
    es <- do.call(rbind, es)
    
    es <- as(es, "dgCMatrix")
    
  } else {
    
    Z <- t(scale(t(X)))
    
    es <- bplapply(geneSets, pcavectorgset, Z,
                   BPPARAM=BPPARAM)
    
    es <- do.call(rbind, es)
    
    if (length(geneSets) == 1)
      es <- matrix(es, nrow=1)
    
    rownames(es) <- names(geneSets)
    colnames(es) <- colnames(X)
  }
  
  es
}


gsva_plage_pca = plage_pca(data.matrix(ds_tmp), gsc)



a <- matrix(rnorm(100000), ncol=20)
pca_out <- BiocSingular::runPCA(a, rank = 30)
pca_1 = pca_out$rotation[ ,1]
pca_1_2 = princomp(a)
pca_new2 = calculatePCA(a, ncomponents = 3,
                        ntop = 500,
                        subset_row = NULL,
                        scale=TRUE,
                        transposed = FALSE)
pca_new_1=pca_new2[,1 ]


svd_out <- BiocSingular::runExactSVD(a)
svd_1 = svd_out$v[, 1]
svd_1_2= svd(a)

umap_out <- umap(a)
tsne = calculateTSNE(a, ncomponents = 3,
                     ntop = 500,
                     subset_row = NULL,
                     scale=TRUE,
                     transposed = FALSE,
                     do.pca = FALSE,
                     seed = 12345)
b = Z[gsc, ]
Z <- Matrix::t(data.matrix(ds_tmp))
Z <- .dgCapply(Z, scale, 2)
Z <- Matrix::t(Z)

umap = calculateUMAP(a, ncomponents = 3,
                     ntop = 500,
                     subset_row = NULL,
                     scale=TRUE,
                     transposed = FALSE)


tsne_out <- tsne(a, max_iter = 10, epoch = 10)



n_cells <- FetchData(X, 
                     vars = c("gene", "gs")) %>%
  group_by(gs) %>%
  dplyr::count(gene) %>% 
  spread(ident, n) 

# View table
View(n_cells)



DimPlot(X,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
X= ScaleData(X)

s_object = CreateSeuratObject(counts = X, min.features  = 10)

spca = RunPCA(object = s_object, features=colnames(X), verbose = FALSE)




s_umap = RunUMAP(s_object, reduction = "pca", dims = 1:20)

VlnPlot(s_object)


DimPlot(object = a, reduction = 'umap')





head(X)







rightsingularsvdvectorgset <- function(gSetIdx, Z) {
  if(is(Z, "dgCMatrix")){
    s <- BiocSingular::runPCA(Z[gSetIdx, ])
  } else {
    s <- BiocSingular::runPCA(Z[gSetIdx, ])
  }
  s$v[, 1]
}


X = data.matrix(ds_tmp)
geneSets = gsc

if(is(X, "dgCMatrix")){
  message("Please bear in mind that this method first scales the values of the gene
    expression data. In order to take advantage of the sparse Matrix type, the scaling
    will only be applied to the non-zero values of the data. This is a provisional 
    solution in order to give support to the dgCMatrix format.")
  
  Z <- Matrix::t(X)
  Z <- .dgCapply(Z, scale, 2)
  Z <- Matrix::t(Z)
  
  es <- lapply(matrix(geneSets), rightsingularsvdvectorgset, Z)
  
  es <- do.call(rbind, es)
  
  es <- as(es, "dgCMatrix")
  
} else {
  
  Z <- t(scale(t(X)))
  
  es <- lapply(matrix(geneSets), rightsingularsvdvectorgset, Z)
  
  es <- do.call(rbind, es)
  
  if (length(geneSets) == 1)
    es <- matrix(es, nrow=1)
  
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
}


#### compare methods


gsva_pval = t.test()





