

install.packages("devtools")
library(devtools)
install_github("rcastelo/GSVA")




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("EnrichmentBrowser")


install.packages("BiocManager")
BiocManager::install("GSVA", version = "devel")

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

BiocManager::install("BiocParallel")


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

library("BiocParallel")
library('data.table')
library("dplyr") 

setwd('C:/Users/zelda/mgr')




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


data<-read.table("Info_files/KEGG_GS.txt",sep='\t')
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
##for (nset in 1:2){
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
## for (dset in 1:2){
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
  ## order_var_tmp = binded_tmp[order(binded_tmp$ls_symbols_tmp, -(as.numeric(binded_tmp$rows_var)[as.numeric(binded_tmp$rows_var)])),]
  order_var_tmp = binded_tmp[with(binded_tmp, order(ls_symbols_tmp,rows_var)),]
  ## remove duplicates
  unique_f_tmp = order_var_tmp[!duplicated(order_var_tmp$ls_symbols_tmp),]
  file_name <- paste("Results/ds_expr_enterez2/ds_uniq", "_",dset,"_e.csv", sep="")
  unique_f_tmp = subset(unique_f_tmp, select=-c(rows_var))
  ## unique_f_tmp = unique_f_tmp[,-2]
  unique_f_tmp = na.omit(unique_f_tmp)
  rownames(unique_f_tmp) = unique_f_tmp$ls_symbols_tmp
  unique_f_tmp[,1] = NULL
  write.csv(unique_f_tmp, file_name, row.names = TRUE)
  ## create object and add to list
  ## object_tmp = new("EachDatasets", dataset_name=colnames(dataset_1_expr)[dset], features_uni=dim(unique_f_tmp), uni_expr=unique_f_tmp)
  
}



## load gene sets 
data<-read.table("Info_files/KEGG_GS.txt",sep='\t')
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

## methods_gsva = c( "zscore", "plage")
methods_gsva = c("zscore")


for (m in methods_gsva){
  
  for (i in 1:length(geo2kegg)){
    # read from file expression set: rows - enterezID, cols - genes
    file_name <- paste("Results/ds_expr_enterez2/ds_uniq", "_",i,"_e.csv", sep="")
    ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
    # main function
    ES_GSVA = gsva(data.matrix(ds_tmp), gsc,
                   method= "zscore",
                   kcdf="Gaussian",
                   abs.ranking=FALSE,
                   min.sz=2,
                   max.sz=Inf, ## all paths even if short
                   parallel.sz=1L,
                   mx.diff=TRUE,
                   ssgsea.norm=TRUE,
                   verbose=TRUE,
                   BPPARAM=SnowParam(2, progressbar = FALSE))
    # write to apropriate folder
    file_out <- paste("Results/new_res/", m, "/ds", "_",i,"_e.csv", sep="")
    write.csv(ES_GSVA, file_out, row.names = TRUE)
  }
  m
}


file_name <- paste("Results/ds_expr_enterez2/ds_uniq_42_e.csv", sep="")
ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
ES_GSVA = gsva(data.matrix(ds_tmp), gsc,
               method="zscore",
               kcdf="Gaussian",
               abs.ranking=FALSE,
               min.sz=2,
               max.sz=Inf, ## all paths even if short
               parallel.sz=1L,
               mx.diff=TRUE,
               ssgsea.norm=TRUE,
               verbose=TRUE)



### calculate t-test


## make df for each method

gsva_ttest <- data.frame(matrix(ncol = 5, nrow = 42))
colnames(gsva_ttest) <- c("DataSet", "t_value","p_value", 'c_nbr', 'd_nbr')

plage_ttest <- data.frame(matrix(ncol = 5, nrow = 42))
colnames(plage_ttest) <- c("DataSet", "t_value","p_value", 'c_nbr', 'd_nbr')

ssgsea_ttest <- data.frame(matrix(ncol = 5, nrow = 42))
colnames(ssgsea_ttest) <- c("DataSet", "t_value","p_value", 'c_nbr', 'd_nbr')

zscore_ttest <- data.frame(matrix(ncol = 5, nrow = 42))
colnames(zscore_ttest) <- c("DataSet", "t_value","p_value", 'c_nbr', 'd_nbr')

path_1 = 'Results/ds_expr_enterez2'

method_ttest <- function(df_ttest, path_f) {
  for (nset in 1:length(geo2kegg)){
  ##for (nset in 1:2){
    ## controles / diseases
    dataset_tmp = geo2kegg[[nset]]
    c_d_data = dataset_tmp@phenoData@data
    nbr_c = count(c_d_data, vars = "Group")
    df_ttest$DataSet[nset] = nset
    df_ttest$d_nbr[nset] = sum(c_d_data$Group=="d")
    df_ttest$c_nbr[nset] = sum(c_d_data$Group=="c")
    
    t_df <- transpose(c_d_data)
    colnames(t_df) <- c_d_data$Sample
    
    file_name <- paste(path_f, 'ds', "_",nset,"_e.csv", sep="")
    ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
    c_sample = c_d_data$Sample[c_d_data$Group== 'c']
    d_sample = c_d_data$Sample[c_d_data$Group== 'd']
    c_mean = mean(as.matrix(ds_tmp[ , c(c_sample)]))
    d_mean = mean(as.matrix(ds_tmp[ , c(d_sample)]))
    
    ## transpose 
    t_tmp = transpose(ds_tmp)
    colnames(t_tmp) <- rownames(ds_tmp)
    rownames(t_tmp) = colnames(ds_tmp)
    
    ## calculate means for each row
    t_tmp$Means = rowMeans(t_tmp[,-1])
    t_tmp$Group = c_d_data$Group
    
    ##  t-test
    df <- t_tmp %>%
      filter(Group == "c" | Group == "d") %>%
      select(Group, Means)
    
    t1 = t.test(Means ~ Group, data = df)
    df_ttest$t_value[nset] = t1$statistic
    df_ttest$p_value[nset] = t1$p.value
    
  }
    
  return(df_ttest)
}



gsva_ttest = method_ttest(gsva_ttest, 'Results/new_res/gsva/')
plage_ttest = method_ttest(plage_ttest, 'Results/new_res/plage/')
ssgsea_ttest = method_ttest(ssgsea_ttest, 'Results/new_res/ssgsea/')

write.csv(gsva_ttest, 'Results/new_res/gsva_ttest.csv', row.names = TRUE)
write.csv(plage_ttest, 'Results/new_res/plage_ttest.csv', row.names = TRUE)
write.csv(ssgsea_ttest, 'Results/new_res/ssgsea_ttest.csv', row.names = TRUE)



df_ttest = gsva_ttest
path_f = 'Results/new_res/gsva/'




gsva_pval <- data.frame(matrix(ncol = 1, nrow = 1))
colnames(gsva_pval) = ('genes')



t_test_all_genes <- function(path_f) {
gsva_pval <- data.frame(matrix(ncol = 1, nrow = 1))
colnames(gsva_pval) = c('genes')

for (nset in 1:length(geo2kegg)){
  ##for (nset in 1:2){
  ## controles / diseases
  dataset_tmp = geo2kegg[[nset]]
  c_d_data = dataset_tmp@phenoData@data
  nbr_c = count(c_d_data, vars = "Group")
  df_ttest$DataSet[nset] = nset
  df_ttest$d_nbr[nset] = sum(c_d_data$Group=="d")
  df_ttest$c_nbr[nset] = sum(c_d_data$Group=="c")
  
  t_df <- transpose(c_d_data)
  colnames(t_df) <- c_d_data$Sample
  
  file_name <- paste(path_f, 'ds', "_",nset,"_e.csv", sep="")
  ds_tmp = read.csv(file_name,check.names=FALSE, row.names=1)
  c_sample = c_d_data$Sample[c_d_data$Group== 'c']
  d_sample = c_d_data$Sample[c_d_data$Group== 'd']
  ds_tmp$mean_c = rowMeans(as.data.frame(ds_tmp[ , c(c_sample)])[,-1])
  ds_tmp$mean_d = rowMeans(as.data.frame(ds_tmp[ , c(d_sample)])[,-1])
  #ds_tmp$p_value = apply(ds_tmp, 1, function(x) t.test(x[-3],x[-2])$p.value)
  ds_tmp$p_value2 = apply(ds_tmp, 1, function(x) t.test(x[c(c_sample)],x[c(d_sample)])$p.value)
  ds_tmp$genes = rownames(ds_tmp)
  p_frame = as.data.frame(ds_tmp[ ,c('genes', 'p_value2') ])
  colnames(p_frame) = c('genes', nset)
  gsva_pval = merge(p_frame, gsva_pval, by= "genes",all = TRUE)
  
  
}
return(gsva_pval)
}



gsva_p = t_test_all_genes('Results/new_res/gsva/')
plage_p = t_test_all_genes('Results/new_res/plage/')
ssgsea_p = t_test_all_genes('Results/new_res/ssgsea/')

write.csv(gsva_p, 'Results/new_res/gsva_all_genes.csv', row.names = TRUE)
write.csv(plage_p, 'Results/new_res/plage_all_genes.csv', row.names = TRUE)
write.csv(ssgsea_p, 'Results/new_res/ssgsea_all_genes.csv', row.names = TRUE)



a = as.data.frame(a)
a$mean_c = rowMeans(a[,-1])
pValues <- apply(ds_tmp, 1, function(x) t.test(x[c(c_sample)],x[c(d_sample)])$p.value)

t.test(ds_tmp$mean_c[1],ds_tmp$mean_d[1])



h1 = ds_tmp[ds_tmp$p_value2 < 0.05, ]
gsva_pval = rbind(gsva_pval,as.data.frame(ds_tmp[ , 'p_value2']))
gsva_pval = merge(p_frame, gsva_pval, by= "genes",all.x = TRUE)
