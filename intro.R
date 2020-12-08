



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("EnrichmentBrowser")

BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133plus2.db")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("KEGGREST")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GSEABenchmarkeR")
BiocManager::install("DataCombine")

library("DataCombine")
library("KEGGREST")
library("GSEABenchmarkeR")
library("EnrichmentBrowser")
library("GSVA")
library("GSEABase")
library("AnnotationDbi")
library("hgu133a.db")    ##for Human
library("hgu133plus2.db")
library("plyr")
library("org.Hs.eg.db")

setwd('C:/Users/hp/Desktop/mgr')




## load data

geo2kegg = loadEData("geo2kegg")
## 42 datasety
## assayData: x features, y samples 
## - hat contains expression levels of x probe sets measured for y patients.


## summarize expression levels for probes annotated to the same gene

geo2kegg_sum = maPreproc(geo2kegg)
names(geo2kegg)


## make df with diseases
df_diseases <- data.frame(matrix(ncol = 6, nrow = 42))
colnames(df_diseases) <- c("ID", "Disease","Features","Samples","Annotation","TargetGenSet")

for (i in 1:length(geo2kegg)){
  dataset_tmp = geo2kegg[[i]]
  df_diseases$ID[i] = dataset_tmp@experimentData@name
  df_diseases$Disease[i] = dataset_tmp@experimentData@other$disease
  df_diseases$Features[i] = dim(dataset_tmp)[1]
  df_diseases$Samples[i] = dim(dataset_tmp)[2]
  df_diseases$Annotation[i] = dataset_tmp@annotation
  df_diseases$TargetGenSet[i] = dataset_1@experimentData@other$targetGeneSets
  
  
  
}

unique_disease = unique(lapply(df_diseases$Disease, tolower))
write.csv(unique_disease,"diseases.csv", row.names = TRUE)
write.csv(df_diseases,"KEGG_diseases.csv", row.names = TRUE)

## get disease code
data.dir <- system.file("extdata", package="GSEABenchmarkeR")
mala.kegg.file <- file.path(data.dir, "malacards", "KEGG.rds")
mala.kegg <- readRDS(mala.kegg.file)
sapply(mala.kegg, nrow)

## Mapping between dataset ID and disease code
d2d.file <- file.path(data.dir, "malacards", "GseId2Disease.txt")
d2d.map <- readDataId2diseaseCodeMap(d2d.file)
head(d2d.map)


res <- keggInfo("kegg") ## displays the current statistics of the KEGG database
cat(res)
res <- keggInfo("pathway") 
## displays the number pathway entries including both
## the reference and organism-specific pathways
cat(res)
res <- keggLink("pathway", "hsa") ## displays the number of gene entries for the
## KEGG organism Homo sapiens
cat(res)


# get info about 1 path:
info_hsa05418 <- keggLink("hsa05418")
get_info_hsa05418 <- keggGet("hsa05418")
get_hsa10000 = keggGet("hsa:10000")
## estimates GSVA enrichment scores:


gs_1 = g1$GSE1297$Sample
geneSets <- list(set1=paste("g", 1:3, sep=""),
                 set2=paste("g", 4:6, sep=""),
                 set3=paste("g", 7:10, sep=""))


## new keegs
listDatabases()
hsa_list = keggList("hsa")


data<-read.table("KEGG_GS.txt",sep='\t')
## GeneSetCollection KEGGREST
Gz<-list()
for (i in 1:dim(data)[1]){
  test<-as.numeric(data[i,3:length(data)])
  test<-as.character(na.omit(test))
  Gz[[i]]<-GeneSet(geneIds=test, setName=as.character(data[i,1]),collectionType=KEGGCollection(as.character(data[i,1])), setIdentifier=as.character(data[i,1]))
}
gsc <- GeneSetCollection(Gz)


## get one expr set
dataset_1=geo2kegg[[1]]
dataset_1_expr = dataset_1@assayData$exprs
write.csv(dataset_1_expr,"dataset_1_expr.csv", row.names = TRUE)
all_rows_expr = rownames(dataset_1_expr)
unique_rows_expr = unique(all_rows_expr)


### get gen symbol for each of probeid
select(org.Hs.eg.db, c("222247_at"), c("SYMBOL","ENTREZID", "GENENAME")) ## trying example

PROBES<- as.character(dataset_1_expr[,1, drop=FALSE])
OUT <- select(hgu133a.db,rownames(dataset_1_expr), columns=c("SYMBOL", "ENTREZID", "GENENAME"))

PROBESplus2<- as.character(dataset_1_expr[,1, drop=FALSE])
OUTplus2 <- select(hgu133plus2.db,rownames(dataset_1_expr), columns=c("SYMBOL", "ENTREZID", "GENENAME"))

write.csv(OUT,"OUT.csv", row.names = TRUE)


## list of symbols
list_symbols = OUT$SYMBOL
ls_symbols = as.list(list_symbols)


## make df with expressions, symbols and var
rows_var = rowVars(dataset_1_expr)
dataset_1_expr_symbols = cbind(dataset_1_expr,rows_var,list_symbols)
dataset_1_expr_symbols = as.data.frame(dataset_1_expr_symbols)
write.csv(dataset_1_expr_symbols,"dataset_1_expr_sym.csv", row.names = TRUE)


## order by var
order_var = dataset_1_expr_symbols[order(dataset_1_expr_symbols$list_symbols, -(as.numeric(levels(dataset_1_expr_symbols$rows_var))[dataset_1_expr_symbols$rows_var])),]
## remove duplicates
unique_expr_symbols = order_var[!duplicated(order_var$list_symbols),]

setClass("EachDatasets", slots=list(dataset_name="character", features_uni="numeric", uni_expr="list"))
ds1_object = new("EachDatasets", dataset_name=colnames(dataset_1_expr)[1], features_uni=dim(unique_expr_symbols), uni_expr=unique_expr_symbols)

## remove 2 last columns
unique_expr_symbols = subset(unique_expr_symbols, select = -c(rows_var,list_symbols) )

## save only unique expresions
write.csv(unique_expr_symbols, "ds_expr_uni/ds_1_uniq.csv", row.names = TRUE)
## ls_ds_objects = list()
## ls_ds_objects[1] = ds1_object


## make df with nymber of controls and diseases
df_c_d <- data.frame(matrix(ncol = 3, nrow = 42))
colnames(df_c_d) <- c("DataSet", "Diseases","Control")

for (nset in 1:length(geo2kegg)){
  ## controles / diseases
  dataset_tmp = geo2kegg[[nset]]
  c_d_data = dataset_tmp@phenoData@data
  nbr_c = count(c_d_data, vars = "Group")
  df_c_d$DataSet[nset] = nset
  df_c_d$Diseases[nset] = nbr_c$freq[nbr_c$Group=="d"]
  df_c_d$Control[nset] = nbr_c$freq[nbr_c$Group=="c"]
  
}
write.csv(df_c_d, "control_diseases.csv", row.names = TRUE)  


### Function to delete NA  
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}


## create dataset expr unique for each dataset --- update save with enterez id as rownames!!!!
for (dset in 1:length(geo2kegg)){
  dataset_tmp=geo2kegg[[dset]]
  dataset_tmp_expr = dataset_tmp@assayData$exprs
  probes_tmp<- as.character(dataset_tmp_expr[,1, drop=FALSE])
  out_tmp = select(hgu133a.db,rownames(dataset_tmp_expr), columns=c("SYMBOL", "ENTREZID", "GENENAME"))
  symbols_tmp = out_tmp$ENTREZID
#  ls_symbols_tmp = as.list(symbols_tmp)
  ls_symbols_tmp = symbols_tmp
  rows_var = rowVars(dataset_tmp_expr)
#  lengths = max(c(length(dataset_tmp_expr), length(symbols_tmp), length(ls_symbols_tmp)))
#  length(dataset_tmp_expr) = lengths
#  length(symbols_tmp) = lengths
#  length(ls_symbols_tmp) = lengths
  binded_tmp = cbind(ls_symbols_tmp, dataset_tmp_expr, rows_var)
  binded_tmp = as.data.frame(binded_tmp)
  
  ## order by var
  order_var_tmp = binded_tmp[order(binded_tmp$ls_symbols_tmp, -(as.numeric(binded_tmp$rows_var)[as.numeric(binded_tmp$rows_var)])),]
  ## remove duplicates
  unique_f_tmp = order_var_tmp[!duplicated(order_var_tmp$ls_symbols_tmp),]
  file_name <- paste("ds_expr_enterez/ds_uniq", "_",dset,".csv", sep="")
  unique_f_tmp = subset(unique_f_tmp, select=-c(rows_var))
  unique_f_tmp = unique_f_tmp[,-2]
  unique_f_tmp = delete.na(unique_f_tmp)
  rownames(unique_f_tmp) = unique_f_tmp$ls_symbols_tmp
  unique_f_tmp[,1] = NULL
  write.csv(unique_f_tmp, file_name, row.names = TRUE)
  ## create object and add to list
  ## object_tmp = new("EachDatasets", dataset_name=colnames(dataset_1_expr)[dset], features_uni=dim(unique_f_tmp), uni_expr=unique_f_tmp)
  
}



## load gene sets 
data<-read.table("KEGG_GS.txt",sep='\t')
## GeneSetCollection KEGGREST
Gz<-list()
for (i in 1:dim(data)[1]){
  test<-as.numeric(data[i,3:length(data)])
  test<-as.character(na.omit(test))
  Gz[[i]]<-GeneSet(geneIds=test, setName=as.character(data[i,1]),collectionType=KEGGCollection(as.character(data[i,1])), setIdentifier=as.character(data[i,1]))
}
gsc <- GeneSetCollection(Gz)



g_nbr = 1
data_sets = list()
for (gg in geo2kegg){
  set_name = paste("set", g_nbr, sep="")
  data_sets[[set_name]] = c( gg$Sample)
  g_nbr = g_nbr+1
}


## load gen sets from csv
genSets_file = read.csv("KEGG_GS_tocsv.csv", fill=TRUE)
for (i in 1:length(genSets_file$has)){
  name = genSets_file$has[i]
  times0 = 5-nchar(as.character((name)))
  new_name = paste("hsa", strrep("0",times0), sep="")
  new_name = paste(new_name, name, sep="")
  genSets_file$has[i] = new_name
}

for (i in 1:length(genSets_file$has)){
  rownames(genSets_file) <- genSets_file[,1]
  genSets_file[,1] = NULL
}



###################################################################
################## R U N     G S V A ##############################
## methods "gsva", "ssgsea", "zscore", "plage"




## get enterez name for each gse
test_OUT <- select(hgu133a.db,rownames(test), columns=c("SYMBOL", "ENTREZID", "GENENAME"))
## remove duplication of enterez 
test_OUT_uni = test_OUT[!duplicated(test_OUT$PROBEID),]
rownames(test) = test_OUT_uni$ENTREZID
row.names(test) = as.character(test_OUT_uni$ENTREZID)
row.names(test) <- gsub("[a-zA-Z ]", "",  row.names(test))
.rowNamesDF(test, make.names=FALSE) <- test_OUT_uni$ENTREZID
test =  data.matrix(test)
write.csv(test, "test.csv", row.names = TRUE)  



### annotation:  "org.Hs.eg.db" 
### run for each method

methods_gsva = c("gsva","ssgsea", "zscore", "plage")

for (m in methods_gsva){
  
  for (i in 1:length(geo2kegg)){
    # read from file expression set: rows - enterezID, cols - genes
    file_name <- paste("ds_expr_enterez/ds_uniq", "_",i,".csv", sep="")
    ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
    # main function
    ES_GSVA = gsva(data.matrix(ds_tmp), gsc,
                   method= m,
                   kcdf="Gaussian",
                   abs.ranking=FALSE,
                   min.sz=1,
                   max.sz=Inf, ## all paths even if short
                   parallel.sz=2,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   verbose=TRUE)
    # write to apropriate folder
    file_out <- paste(m, "/ds", "_",i,".csv", sep="")
    write.csv(ES_GSVA, file_out, row.names = TRUE)
  }
  m
}



ES_GSVA = gsva(data.matrix(ds1_uni), gsc,
               method="gsva",
               kcdf="Gaussian",
               abs.ranking=FALSE,
               min.sz=1,
               max.sz=Inf, ## all paths even if short
               parallel.sz=1L,
               mx.diff=TRUE,
               ssgsea.norm=TRUE,
               verbose=TRUE)


