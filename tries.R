

BiocManager::install(version = '4.0')

BiocManager::install("GSVA", version = "devel")
BiocManager::install("Seurat")

install.packages("remotes")
library(remotes)
install_github("salmonellix/GSVA", force=TRUE)
library("GSVA")
library("BiocParallel")
library("BiocSingular")

remotes::install_github("davismcc/scater")
library("scater")
library("davismcc/scater")

BiocManager::install("scater")
install.packages("umap")
library("umap")

install.packages("tsne")
library("tsne")
library("dplyr")
library("tidyr")
library("Seurat")

methods_gsva = c("zscore")








for (m in methods_gsva){
  
  for (i in 1:length(geo2kegg)){
    # read from file expression set: rows - enterezID, cols - genes
    file_name <- paste("Results/ds_expr_enterez2/ds_uniq", "_",i,"_e.csv", sep="")
    ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
    # main function
    ES_GSVA = gsva(data.matrix(ds_tmp), gsc,
                   method= c("zscore"),
                   kcdf=c("Gaussian"),
                   abs.ranking=FALSE,
                   min.sz=1,
                   max.sz=Inf, ## all paths even if short
                   parallel.sz=1L,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   verbose=FALSE)
    # write to apropriate folder
    file_out <- paste("Results/modifications/zscore_fisher//ds", "_",i,"_zf.csv", sep="")
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
pca_new2 = calculatePCA(data.matrix(a), ncomponents = 20,
                        ntop = 500,
                        subset_row = NULL,
                        scale=TRUE,
                        transposed = FALSE)
pca_new_1=pca_new2[,1 ]


svd_out <- BiocSingular::runExactSVD(a)
svd_1 = svd_out$v[, 1]
svd_1_2= svd(a)

umap_out <- umap(a)
tsne = calculateTSNE(a)
umap = calculateUMAP(a)


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

es