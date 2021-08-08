library(ggplot2)
library(ClustImpute)
library(ClusterR)
library(mice)
library(Amelia)
library(missRanger)
library(microbenchmark)
library(cowplot)
library(dplyr)
library(tidyr)
library(Hmisc)
library(corrplot)

data("iris")
dat <- as.data.frame(scale(iris[,1:4]))
true_label <- as.numeric(iris$Species)

type_missing <- "MAR" # use MCAR or MAR

### Parameters
# ClustImpute
nr_iter <- 10 # iterations of procedure
n_end <- 6 # step until convergence of weight function to 1; set to 1 to eliminate the effect of the weight function
nr_cluster <- 3 # number of clusters
c_steps <- 10 # numer of cluster steps per iteration
# Clustering based on imputation
# nr_iter_other <- (nr_iter-n_end) * c_steps # comparable number of steps "after" imputation
nr_iter_other <- nr_iter * c_steps

param_times <- 30 # how often to compute the benchmark

### for various p
P <- c(.05,.1,.2,.3,.4,.5,.6,.7)

results <- list()
count <- 0

set.seed(124) # for reproducibility
for (p in P) {
  count <- count + 1
  # create missings for data
  dat_with_miss <- miss_sim(dat,p,seed_nr=739,type=type_missing)
  
  png(paste0("Corrplot missings ",type_missing," ",p,".png"))
  mis_ind <- is.na(dat_with_miss)
  corrplot(cor(mis_ind),method="number")
  dev.off()
  
  ### Benchmark time and calculate rand index
  rand_ClustImpute <- 0
  rand_ClustImpute_false <- 0
  rand_ClustImpute_no_weight <- 0
  rand_missRanger <- 0
  rand_RandomImp <- 0
  rand_mice_pmm <- 0
  rand_mice_cart <- 0
  rand_amelia <- 0
  
  # Use a different seed each time
  mbm <- microbenchmark("ClustImpute" = {
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end, seed_nr = random_seed)
    rand_ClustImpute <- rand_ClustImpute + external_validation(true_label, res$clusters)
  },"ClustImpute_assign_with_wf_false" = {
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end, seed_nr = random_seed, 
                       assign_with_wf =  FALSE)
    rand_ClustImpute_false <- rand_ClustImpute_false + external_validation(true_label, res$clusters)
  },
  "ClustImpute_no_weight" = {
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=1, seed_nr = random_seed)
    rand_ClustImpute_no_weight <- rand_ClustImpute_no_weight + external_validation(true_label, res$clusters)
  },
  "RandomImp"={
    dat_random_imp <- dat_with_miss
    for (j in 1:dim(dat)[2]) {
      dat_random_imp[,j] <- impute(dat_random_imp[,j],fun="random")
    }
    cl_RandomImp <- KMeans_arma(data=dat_random_imp,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
    pred_RandomImp <- predict_KMeans(dat_random_imp,cl_RandomImp)
    class(pred_RandomImp) <- "numeric"
    rand_RandomImp <- rand_RandomImp + external_validation(true_label, pred_RandomImp)
  },
  "MissRanger"={
    imp_missRanger <- missRanger(dat_with_miss,pmm.k = 3,verbose = 0)
    cl_missRanger <- KMeans_arma(data=imp_missRanger,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
    pred_missRanger <- predict_KMeans(imp_missRanger,cl_missRanger)
    class(pred_missRanger) <- "numeric"
    rand_missRanger <- rand_missRanger + external_validation(true_label, pred_missRanger)
  }, 
  "MICE_pmm"={
    dat_mice <- mice(dat_with_miss,m=1,maxit=50,meth='pmm',printFlag=FALSE) # single data set
    dat_mice <- mice::complete(dat_mice)
    if (sum(is.na(dat_mice))==0) {
      cl_mice <- KMeans_arma(data=dat_mice,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
      pred_mice <- predict_KMeans(dat_mice,cl_mice)
      class(pred_mice) <- "numeric"
      rand_mice_pmm <- rand_mice_pmm + external_validation(true_label, pred_mice)
    } else
      rand_mice_pmm <- NA
  },
  "MICE_cart"={
    dat_mice <- mice(dat_with_miss,m=1,maxit=50,meth='cart',printFlag=FALSE) # single data set
    dat_mice <- mice::complete(dat_mice)
    if (sum(is.na(dat_mice))==0) {
      cl_mice <- KMeans_arma(data=dat_mice,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
      pred_mice <- predict_KMeans(dat_mice,cl_mice)
      class(pred_mice) <- "numeric"
      rand_mice_cart <- rand_mice_cart + external_validation(true_label, pred_mice)
    } else
      rand_mice_cart <- NA
  }, 
  "AMELIA"= {
    if (p<.7) {
      dat_amelia <- amelia(dat_with_miss,m=1,p2s=0) # single data set
      dat_amelia <- dat_amelia$imputations$imp1
      if (sum(is.na(dat_amelia))==0) {
        cl_amelia <- KMeans_arma(data=dat_amelia,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
        pred_amelia <- predict_KMeans(dat_amelia,cl_amelia)
        class(pred_amelia) <- "numeric"
        rand_amelia <- rand_amelia + external_validation(true_label, pred_amelia)
      } else
        rand_amelia = NA
    } else rand_amelia <- NA
  }  ,times = param_times,setup = {random_seed=round(runif(1)*1e5)}, unit = "s") # random seed improves performance of ClustImpute
  
  # compute average rand index
  rand_ClustImpute <- rand_ClustImpute/param_times
  rand_ClustImpute_false <- rand_ClustImpute_false/param_times
  rand_ClustImpute_no_weight <- rand_ClustImpute_no_weight/param_times
  rand_missRanger <- rand_missRanger/param_times
  rand_RandomImp <- rand_RandomImp/param_times
  rand_mice_pmm <- rand_mice_pmm/param_times
  rand_mice_cart <- rand_mice_cart/param_times
  rand_amelia <- rand_amelia/param_times
  
  results$randIndex[[count]] <- c(ClustImpute=rand_ClustImpute,
                                  ClustImpute_assign_with_wf_false = rand_ClustImpute_false,
                                  ClustImpute_no_weight=rand_ClustImpute_no_weight,
                                  RandomImp=rand_RandomImp,
                                  missRanger=rand_missRanger,
                                  MICE_pmm=rand_mice_pmm,
                                  MICE_cart=rand_mice_cart,
                                  AMELIA=rand_amelia)

  results$benchmark[[count]] <- mbm
  
  mbm_median <- print(mbm)$median
  results$benchmark_median[[count]] <- mbm_median
}

### Rand index
randtbl <- data.frame(matrix(unlist(results$randIndex),nrow=length(results$randIndex), byrow=T))
colnames(randtbl) <- names(results$randIndex[[1]])
randtbl$p <- P
randtbl

write.csv(randtbl,paste(type_missing,"-","Rand indices iris.csv"))


### plot for rand index
data2plot <- randtbl %>% pivot_longer(-p,names_to = "Algorithm",values_to = "Rand_Index")
rand_plot <- ggplot(data2plot,aes(x=p,y=Rand_Index,color=Algorithm)) + geom_point() + theme_cowplot() + geom_line() +
  xlab("Share of missings") + ylab("Rand Index") + scale_x_continuous(breaks = 1:7/10) + scale_y_continuous(breaks = 0:10/10) + 
  theme(legend.position="bottom")
rand_plot
save_plot(paste(type_missing,"-","Rand index iris.png"),rand_plot,base_aspect_ratio = 2.5)


### Violinplot performance
violin_plot <- plot_grid(autoplot(results$benchmark[[1]]),autoplot(results$benchmark[[2]]),autoplot(results$benchmark[[3]]),
                         autoplot(results$benchmark[[4]]),autoplot(results$benchmark[[5]]),autoplot(results$benchmark[[6]]),
                         autoplot(results$benchmark[[7]]),autoplot(results$benchmark[[8]]),autoplot(results$benchmark[[9]]),
               ncol=2,labels=sprintf("p = %s",P),label_size = 12,vjust=1)
violin_plot
save_plot(paste(type_missing,"-","Violinplot running time iris.png"),violin_plot,base_aspect_ratio = 2.5)


# save results
save(results,rand_plot,violin_plot,file=paste(type_missing,"-","results_benchmarking_iris.Rdata"))

