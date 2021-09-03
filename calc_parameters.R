
library("ggplot2")
library("scater")
library(reshape2)

setwd('C:/Users/zelda/mgr')




runSim <- function(p, n, gs.sz, S2N, fracDEgs) {
  sizeDEgs <- round(fracDEgs * gs.sz)
  group.n <- round(n / 2)
  sampleEffect <- rnorm(n, mean=0, sd=1)
  sampleEffectDE <- rnorm(n, mean=S2N, sd=0.5)
  probeEffect <- rnorm(p, mean=0, sd=1)
  noise <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  noiseDE <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  M <- outer(probeEffect, sampleEffect, "+") + noise
  M2 <- outer(probeEffect, sampleEffectDE, "+") + noiseDE
  M[1:sizeDEgs, 1:group.n] <- M2[1:sizeDEgs, 1:group.n]
  rownames(M) <- paste0("g", 1:nrow(M))
  geneSets <- list(H1GeneSet=paste0("g", 1:(gs.sz)),H0GeneSet=paste0("g", (gs.sz+1):(2*gs.sz)))

  es.gsva <- gsva(M, geneSets, method='gsva', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.ss <- gsva(M, geneSets, method='ssgsea',abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.z <- gsva(M, geneSets,method='zscore', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.z_fisher <- gsva(M, geneSets,method='zscore_fisher', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.z_stouffer <- gsva(M, geneSets,method='zscore_stouffer', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.plage <- gsva(M, geneSets,method='plage', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.plage_pca <- gsva(M, geneSets,method='plage_pca', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.plage_umap <- gsva(M, geneSets,method='plage_umap', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  es.plage_tsne <- gsva(M, geneSets,method='plage_tsne',abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  
  h1.gsva.pval <- t.test(es.gsva["H1GeneSet", 1:group.n],es.gsva["H1GeneSet", (group.n+1):n])$p.value
  h1.ssgsea.pval <- t.test(es.ss["H1GeneSet", 1:group.n],es.ss["H1GeneSet", (group.n+1):n])$p.value
  h1.zscore.pval <- t.test(es.z["H1GeneSet", 1:group.n],es.z["H1GeneSet", (group.n+1):n])$p.value
  h1.zscore_fisher.pval <- t.test(es.z_fisher["H1GeneSet", 1:group.n],es.z_fisher["H1GeneSet", (group.n+1):n])$p.value
  h1.zscore_stouffer.pval <- t.test(es.z_stouffer["H1GeneSet", 1:group.n],es.z_stouffer["H1GeneSet", (group.n+1):n])$p.value
  h1.plage.pval <- t.test(es.plage["H1GeneSet", 1:group.n],es.plage["H1GeneSet", (group.n+1):n])$p.value
  h1.plage_pca.pval <- t.test(es.plage_pca["H1GeneSet", 1:group.n],es.plage_pca["H1GeneSet", (group.n+1):n])$p.value
  h1.plage_umap.pval <- t.test(es.plage_umap["H1GeneSet", 1:group.n],es.plage_umap["H1GeneSet", (group.n+1):n])$p.value
  h1.plage_tsne.pval <- t.test(es.plage_tsne["H1GeneSet", 1:group.n],es.plage_tsne["H1GeneSet", (group.n+1):n])$p.value
  
  
  h0.gsva.pval <- t.test(es.gsva["H0GeneSet", 1:group.n],es.gsva["H0GeneSet", (group.n+1):n])$p.value
  h0.ssgsea.pval <- t.test(es.ss["H0GeneSet", 1:group.n],es.ss["H0GeneSet", (group.n+1):n])$p.value
  h0.zscore.pval <- t.test(es.z["H0GeneSet", 1:group.n],es.z["H0GeneSet", (group.n+1):n])$p.value
  h0.plage.pval <- t.test(es.plage["H0GeneSet", 1:group.n],es.plage["H0GeneSet", (group.n+1):n])$p.value
  h0.zscore_fisher.pval <- t.test(es.z_fisher["H0GeneSet", 1:group.n],es.z_fisher["H0GeneSet", (group.n+1):n])$p.value
  h0.zscore_stouffer.pval <- t.test(es.z_stouffer["H0GeneSet", 1:group.n],es.z_stouffer["H0GeneSet", (group.n+1):n])$p.value
  h0.plage_pca.pval <- t.test(es.plage_pca["H0GeneSet", 1:group.n],es.plage_pca["H0GeneSet", (group.n+1):n])$p.value
  h0.plage_umap.pval <- t.test(es.plage_umap["H0GeneSet", 1:group.n],es.plage_umap["H0GeneSet", (group.n+1):n])$p.value
  h0.plage_tsne.pval <- t.test(es.plage_tsne["H0GeneSet", 1:group.n],es.plage_tsne["H0GeneSet", (group.n+1):n])$p.value
  
  c(h1.gsva.pval, h1.ssgsea.pval, h1.zscore.pval, h1.plage.pval,h1.zscore_fisher.pval,h1.zscore_stouffer.pval,h1.plage_pca.pval, h1.plage_umap.pval, h1.plage_tsne.pval,
    h0.gsva.pval, h0.ssgsea.pval, h0.zscore.pval, h0.plage.pval,h0.zscore_fisher.pval,h0.zscore_stouffer.pval,h0.plage_pca.pval, h0.plage_umap.pval, h0.plage_tsne.pval)
}

estPwrTypIerr <- function(pvals, alpha=0.05) {
  N <- ncol(pvals)
  c(1 - sum(pvals[1, ] > alpha)/N, 1 - sum(pvals[2, ] > alpha)/N,1 - sum(pvals[3, ] > alpha)/N, 1 - sum(pvals[4, ] > alpha)/N, 1 - sum(pvals[5, ] > alpha)/N, 1 - sum(pvals[6, ] > alpha)/N,1 - sum(pvals[7, ] > alpha)/N, 1 - sum(pvals[8, ] > alpha)/N,1 - sum(pvals[9, ] > alpha)/N,
    sum(pvals[10, ] <= alpha)/N, sum(pvals[11, ] <= alpha)/N, sum(pvals[12, ] <= alpha)/N, sum(pvals[13, ] <= alpha)/N, sum(pvals[14, ] <= alpha)/N, sum(pvals[15, ] <= alpha)/N, sum(pvals[16, ] <= alpha)/N, sum(pvals[17, ] <= alpha)/N, sum(pvals[18, ] <= alpha)/N)
  }


set.seed(1234)
exp1 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
                + estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
                + estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
                + estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.5))))

exp2 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
                + estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
                + estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
                + estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.5))))
exp3 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
                + estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
                + estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
                + estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.8))))
exp4 <- cbind(estPwrTypIerr(replicate(60, runSim(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
                + estPwrTypIerr(replicate(60, runSim(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
                + estPwrTypIerr(replicate(60, runSim(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
                + estPwrTypIerr(replicate(60, runSim(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.8))))


methods_name = c('gsva', 'ssgsea', 'zscore', 'zscore_fisher', 'zscore_stouffer', 'plage', 'plage_pca', 'plage_umap', 
            'plage_tsne','gsva', 'ssgsea', 'zscore', 'zscore_fisher', 'zscore_stouffer', 'plage', 'plage_pca', 'plage_umap', 'plage_tsne')

exp1=cbind(exp1,methods_name)
exp2=cbind(exp2,methods_name)
exp3=cbind(exp3,methods_name)
exp4=cbind(exp4,methods_name)


colnames(exp4)[1]="n_10"
colnames(exp4)[2]="n_20"
colnames(exp4)[3]="n_40"
colnames(exp4)[4]="n_60"
exp2_stat_power$n_60 <-as.numeric(exp2_stat_power$n_60)
exp2_stat_power$methods_name <-as.character(exp2_stat_power$methods_name)

exp2_stat_power$n_10=round(exp1_stat_power$n_10,3)

exp1_stat_power=exp1[c(1:9),]
exp1_stat_power = as.data.frame(exp1_stat_power)
exp2_stat_power=exp2[c(1:9),]
exp2_stat_power = as.data.frame(exp2_stat_power)
exp3_stat_power=exp3[c(1:9),]
exp3_stat_power = as.data.frame(exp3_stat_power)
exp4_stat_power=exp4
exp4_stat_power = as.data.frame(exp4_stat_power)


exp4_stat_power$n_60 <-as.numeric(exp4_stat_power$n_60)
exp4_stat_power$methods_name <-as.character(exp4_stat_power$methods_name)

exp4_stat_power$n_10=round(exp4_stat_power$n_10,3)

df <- melt(exp4_stat_power, id.vars = "methods_name")  #the function melt reshapes it from wide to long
df$rowid <- 1:9  #add a rowid identifying variable

df$value <-as.numeric(df$value)
df$value <-round(df$value,3) #the "-1" excludes column 1



ggplot(df, aes(variable, value, group=factor(methods_name))) + geom_line(aes(color=factor(methods_name)), linetype = "dashed")+ geom_point(aes(color=factor(methods_name)))+ ylab("Sample size") + 
  ylab("Empirical Type???I Errorr")+ ggtitle("Experiment 4 -  type-I error rate") + theme(plot.title = element_text(size = 20, face = "bold"))+theme(text = element_text(size = 15))


exp1_err1=exp1[c(10:18),]
exp2_err1=exp2[c(10:18),]
exp3_err1=exp3[c(10:18),]
exp4_err1=exp4[c(10:18),]

exp1_err1 = as.data.frame(exp1_err1)
exp2_err1 = as.data.frame(exp2_err1)
exp3_err1 = as.data.frame(exp3_err1)
exp4_err1 = as.data.frame(exp4_err1)


exp4_err1$n_10 <-as.numeric(exp4_err1$n_10)
exp4_err1$n_20 <-as.numeric(exp4_err1$n_20)
exp4_err1$n_40 <-as.numeric(exp4_err1$n_40)
exp4_err1$n_60 <-as.numeric(exp4_err1$n_60)

exp4_err1$methods_name <-as.character(exp4_err1$methods_name)


df <- melt(exp4_err1, id.vars = "methods_name")  #the function melt reshapes it from wide to long
df$rowid <- 1:9  #add a rowid identifying variable

df$value <-as.numeric(df$value)
df$value <-round(df$value,3) #the "-1" excludes column 1


## plots


write.csv(exp1,"exp1.csv", row.names = TRUE)  
write.csv(exp2,"exp2.csv", row.names = TRUE)  
write.csv(exp3,"exp3.csv", row.names = TRUE)  
write.csv(exp4,"exp4.csv", row.names = TRUE)  

aaa = replicate(60, runSim(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.5))




######## calculate time #####################


runSimTime <- function(p, n, gs.sz, S2N, fracDEgs) {
  sizeDEgs <- round(fracDEgs * gs.sz)
  group.n <- round(n / 2)
  sampleEffect <- rnorm(n, mean=0, sd=1)
  sampleEffectDE <- rnorm(n, mean=S2N, sd=0.5)
  probeEffect <- rnorm(p, mean=0, sd=1)
  noise <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  noiseDE <- matrix(rnorm(p*n, mean=0, sd=1), nrow=p, ncol=n)
  M <- outer(probeEffect, sampleEffect, "+") + noise
  M2 <- outer(probeEffect, sampleEffectDE, "+") + noiseDE
  M[1:sizeDEgs, 1:group.n] <- M2[1:sizeDEgs, 1:group.n]
  rownames(M) <- paste0("g", 1:nrow(M))
  geneSets <- list(H1GeneSet=paste0("g", 1:(gs.sz)),H0GeneSet=paste0("g", (gs.sz+1):(2*gs.sz)))
  
  start_time <- Sys.time()
  es.gsva <- gsva(M, geneSets, method='gsva', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  gsva_time = end_time - start_time
  
  start_time <- Sys.time()
  es.ss <- gsva(M, geneSets, method='ssgsea',abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  ss_time = end_time - start_time
  
  
  start_time <- Sys.time()
  es.z <- gsva(M, geneSets,method='zscore', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  z_time = end_time - start_time
  
  start_time <- Sys.time()
  es.z_fisher <- gsva(M, geneSets,method='zscore_fisher', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  z_fisher_time = end_time - start_time
  
  start_time <- Sys.time()
  es.z_stouffer <- gsva(M, geneSets,method='zscore_stouffer', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  z_stouffer_time = end_time - start_time
  
  start_time <- Sys.time()
  es.plage <- gsva(M, geneSets,method='plage', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  plage_time = end_time - start_time
  
  start_time <- Sys.time()
  es.plage_pca <- gsva(M, geneSets,method='plage_pca', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  plage_pca_time = end_time - start_time
  
  
  start_time_u <- Sys.time()
  es.plage_umap <- gsva(M, geneSets,method='plage_umap', abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time_u<- Sys.time()
  plage_umap_time = end_time_u - start_time_u
  
  
  start_time <- Sys.time()
  es.plage_tsne <- gsva(M, geneSets,method='plage_tsne',abs.ranking=FALSE,min.sz=2,max.sz=Inf,parallel.sz=1L,mx.diff=TRUE,ssgsea.norm=TRUE,verbose=FALSE)
  end_time <- Sys.time()
  plage_tsne_time = end_time - start_time

  
  c(gsva_time, ss_time, z_time,z_fisher_time, z_stouffer_time , plage_time ,plage_pca_time, plage_umap_time, plage_tsne_time)
}



set.seed(1234)
time_e1 <- cbind(replicate(60, runSimTime(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.5)),
              + (replicate(60, runSimTime(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
              + (replicate(60,runSimTime(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.5))),
              + (replicate(60, runSimTime(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.5))))

time_e2 <- cbind((replicate(60, runSimTime(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
              + (replicate(60, runSimTime(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
              + (replicate(60, runSimTime(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.5))),
              + (replicate(60, runSimTime(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.5))))
time_e3 <- cbind((replicate(60,runSimTime(1000, 10, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
              + (replicate(60, runSimTime(1000, 20, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
              + (replicate(60, runSimTime(1000, 40, gs.sz=30, S2N=0.5, fracDEgs=0.8))),
              + (replicate(60, runSimTime(1000, 60, gs.sz=30, S2N=0.5, fracDEgs=0.8))))
time_e4 <- cbind((replicate(60, runSimTime(1000, 10, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
              + (replicate(60, runSimTime(1000, 20, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
              + (replicate(60, runSimTime(1000, 40, gs.sz=30, S2N=1.0, fracDEgs=0.8))),
              + (replicate(60, runSimTime(1000, 60, gs.sz=30, S2N=1.0, fracDEgs=0.8))))


write.csv(time_e1,"Results/time_sim/time_e1.csv", row.names = TRUE)  
write.csv(time_e2,"Results/time_sim/time_e2.csv", row.names = TRUE)  
write.csv(time_e3,"Results/time_sim/time_e3.csv", row.names = TRUE)  
write.csv(time_e4,"Results/time_sim/time_e4.csv", row.names = TRUE) 


### load from files ###

time_e1 = read.table("Results/time_sim/time_e1.csv", header = TRUE, sep=',')
time_e1=as.data.frame(time_e1)
time_e1 = time_e1[-c(1)]
time_e2 = read.table("Results/time_sim/time_e2.csv", header = TRUE, sep=',')
time_e2 = time_e2[-c(1)]
time_e3 = read.table("Results/time_sim/time_e3.csv", header = TRUE, sep=',')
time_e3 = time_e3[-c(1)]
time_e4 = read.table("Results/time_sim/time_e4.csv", header = TRUE, sep=',')
time_e4 = time_e4[-c(1)]


time_table1 = data.frame(matrix(ncol = 8, nrow = 9))
colnames(time_table)=c("n_10_mean", "n_10_sd", "n_20_mean", "n_20_sd", "n_40_mean", "n_40_sd", "n_60_mean", "n_60_sd")

methods_name = c('gsva', 'ssgsea', 'zscore', 'zscore_fisher', 'zscore_stouffer', 'plage', 'plage_pca', 'plage_umap', 
                 'plage_tsne')


create.time.table <- function(time_e){
  time_table1 = data.frame(matrix(ncol = 8, nrow = 9))
  
  for (i in 1:9){
    time_table1$n_10_mean[i] = mean(as.numeric(time_e[i,1:60]))
    time_table1$n_20_mean[i] = mean(as.numeric(time_e[i,61:120]))
    time_table1$n_40_mean[i] = mean(as.numeric(time_e[i,121:180]))
    time_table1$n_60_mean[i] = mean(as.numeric(time_e[i,181:240]))
    
    time_table1$n_10_sd[i] = sd(as.numeric(time_e[i,1:60]))
    time_table1$n_20_sd[i] = sd(as.numeric(time_e[i,61:120]))
    time_table1$n_40_sd[i] = sd(as.numeric(time_e[i,121:180]))
    time_table1$n_60_sd[i] = sd(as.numeric(time_e[i,181:240]))
    
  }
  
  time_table1= time_table1[-c(1:8)]
  
}

t1 = create.time.table(time_e1)
t2 = create.time.table(time_e2)
t3 = create.time.table(time_e3)
t4 = create.time.table(time_e4)


write.csv(t1,"Results/time_sim/time_1_stats.csv", row.names = TRUE)  
write.csv(t2,"Results/time_sim/time_2_stats.csv", row.names = TRUE)  
write.csv(t3,"Results/time_sim/time_3_stats.csv", row.names = TRUE)  
write.csv(t4,"Results/time_sim/time_4_stats.csv", row.names = TRUE) 


qplot(methods_name,as.numeric(t1$n_10_mean))+geom_errorbar(aes(x=methods_name, ymin=t1$n_10_mean-t1$n_10_sd, ymax=t1$n_10_mean+t1$n_10_sd), width=0.25)+ylab("Time [s]") + 
  xlab("Method")+ ggtitle("Performance time for each method n = 10") + theme(plot.title = element_text(size = 20, face = "bold"))+theme(text = element_text(size = 15))


qplot(methods_name,as.numeric(t1$n_20_mean))+geom_errorbar(aes(x=methods_name, ymin=t1$n_20_mean-t1$n_20_sd, ymax=t1$n_20_mean+t1$n_20_sd), width=0.25)+ylab("Time [s]") + 
  xlab("Method")+ ggtitle("Performance time for each method n = 20") + theme(plot.title = element_text(size = 20, face = "bold"))+theme(text = element_text(size = 15))

qplot(methods_name,as.numeric(t1$n_40_mean))+geom_errorbar(aes(x=methods_name, ymin=t1$n_40_mean-t1$n_40_sd, ymax=t1$n_40_mean+t1$n_40_sd), width=0.25)+ylab("Time [s]") + 
  xlab("Method")+ ggtitle("Performance time for each method n = 40") + theme(plot.title = element_text(size = 20, face = "bold"))+theme(text = element_text(size = 15))

qplot(methods_name,as.numeric(t1$n_60_mean))+geom_errorbar(aes(x=methods_name, ymin=t1$n_60_mean-t1$n_60_sd, ymax=t1$n_60_mean+t1$n_60_sd), width=0.25)+ylab("Time [s]") + 
  xlab("Method")+ ggtitle("Performance time for each method n = 60") + theme(plot.title = element_text(size = 20, face = "bold"))+theme(text = element_text(size = 15))

t1$method = methods_name




single_run = runSimTime(100, 10, gs.sz=30, S2N=0.5, fracDEgs=0.5)

t1_means = t1[-c(5,6,7,8)]
colnames(t1_means)[1]="n_10"
colnames(t1_means)[2]="n_20"
colnames(t1_means)[3]="n_40"
colnames(t1_means)[4]="n_60"
df <- melt(t1_means, id.vars = "method")  #the function melt reshapes it from wide to long
df$rowid <- 1:9  #add a rowid identifying variable

df$value <-as.numeric(df$value)
df$value <-round(df$value,3) #the "-1" excludes column 1


ggplot(df, aes(variable, value, group=factor(method))) + geom_line(aes(color=factor(method)), linetype = "dashed")+ geom_point(aes(color=factor(method)))+ xlab("Sample size") + 
  ylab("Mean Time [s]")+ ggtitle("Mean time for different sample size") + theme(plot.title = element_text(size = 20, face = "bold"))+theme(text = element_text(size = 15))


