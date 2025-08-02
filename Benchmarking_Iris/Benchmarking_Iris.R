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
dat <- as.data.frame(scale(iris[,1:4])) # assume scaling is known
true_label <- as.numeric(iris$Species)

type_missing <- "MCAR" # use MCAR or MAR

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
P <- c(.05,.1,.2,.3,.4,.5,.7)

count <- 0
results_mean_randIndex <- list()
results_sd_randIndex <- list()

for (p in P) {
  
  count <- count + 1
  # create missings for data
  dat_with_miss <- miss_sim(dat,p,type=type_missing, seed_nr = 555+round(p*1000))
  
  png(paste0("Corrplot missings ",type_missing," ",p,".png"))
  mis_ind <- is.na(dat_with_miss)
  corrplot(cor(mis_ind),method="number")
  dev.off()
  
  ### Benchmark time and calculate rand index
  rand_ClustImpute <- c()
  rand_ClustImpute_false <- c()
  rand_ClustImpute_no_weight <- c()
  rand_missRanger <- c()
  rand_RandomImp <- c()
  rand_mice_pmm <- c()
  rand_mice_cart <- c()
  rand_amelia <- c()

  # Use a different seed each time
  for (k in 1:param_times) {
    
    random_seed = 44444+round(p*1000)+k
    
    # ClustImpute
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end, seed_nr = random_seed)
    rand_ClustImpute <- c(rand_ClustImpute, external_validation(true_label, res$clusters))
    
    # ClustImpute: weight function is only applied in the centroid computation, but ignored in the cluster assignment.
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end, seed_nr = random_seed, 
                       assign_with_wf =  FALSE)
    rand_ClustImpute_false <- c(rand_ClustImpute_false, external_validation(true_label, res$clusters))
    
    # ClustImpute: No weight function
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=1, seed_nr = random_seed)
    rand_ClustImpute_no_weight <- c(rand_ClustImpute_no_weight, external_validation(true_label, res$clusters))
    
    # Random Imputation
    dat_random_imp <- dat_with_miss
    for (j in 1:dim(dat)[2]) {
      dat_random_imp[,j] <- impute(dat_random_imp[,j],fun="random")
    }
    cl_RandomImp <- KMeans_arma(data=dat_random_imp,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
    pred_RandomImp <- predict_KMeans(dat_random_imp,cl_RandomImp)
    class(pred_RandomImp) <- "numeric"
    rand_RandomImp <- c(rand_RandomImp, external_validation(true_label, pred_RandomImp))
    
    # MissRanger
    imp_missRanger <- missRanger(dat_with_miss,pmm.k = 3,verbose = 0)
    cl_missRanger <- KMeans_arma(data=imp_missRanger,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
    pred_missRanger <- predict_KMeans(imp_missRanger,cl_missRanger)
    class(pred_missRanger) <- "numeric"
    rand_missRanger <- c(rand_missRanger, external_validation(true_label, pred_missRanger))
    
    # MICE PMM
    dat_mice <- mice(dat_with_miss,m=1,maxit=50,meth='pmm',printFlag=FALSE) # single data set
    dat_mice <- mice::complete(dat_mice)
    if (sum(is.na(dat_mice))==0) {
      cl_mice <- KMeans_arma(data=dat_mice,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
      pred_mice <- predict_KMeans(dat_mice,cl_mice)
      class(pred_mice) <- "numeric"
      rand_mice_pmm <- c(rand_mice_pmm, external_validation(true_label, pred_mice))
    } else
      rand_mice_pmm <- c(rand_mice_pmm, NA)
    
    # MICE CART
    dat_mice <- mice(dat_with_miss,m=1,maxit=50,meth='cart',printFlag=FALSE) # single data set
    dat_mice <- mice::complete(dat_mice)
    if (sum(is.na(dat_mice))==0) {
      cl_mice <- KMeans_arma(data=dat_mice,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
      pred_mice <- predict_KMeans(dat_mice,cl_mice)
      class(pred_mice) <- "numeric"
      rand_mice_cart <- c(rand_mice_cart, external_validation(true_label, pred_mice))
    } else
      rand_mice_cart <- c(rand_mice_cart, NA)
    
    # AMELIA
    dat_amelia <- amelia(dat_with_miss,m=1,p2s=0,empri = p/7 * nrow(dat_with_miss)) # from doc: "a reasonable upper bound is around 10 percent of the rows of the data."
    dat_amelia <- dat_amelia$imputations$imp1
    if (sum(is.na(dat_amelia))==0) {
      cl_amelia <- KMeans_arma(data=dat_amelia,clusters=nr_cluster,n_iter=nr_iter_other,seed=random_seed)
      pred_amelia <- predict_KMeans(dat_amelia,cl_amelia)
      class(pred_amelia) <- "numeric"
      rand_amelia <- c(rand_amelia, external_validation(true_label, pred_amelia))
    } else
      rand_amelia = c(rand_amelia, NA)
    
  }
  
  # compute average rand index
  na_rm_setting = TRUE
  
  mean_rand_ClustImpute <- mean(rand_ClustImpute,na.rm = na_rm_setting)
  mean_rand_ClustImpute_false <- mean(rand_ClustImpute_false,na.rm = na_rm_setting)
  mean_rand_ClustImpute_no_weight <- mean(rand_ClustImpute_no_weight,na.rm = na_rm_setting)
  mean_rand_missRanger <- mean(rand_missRanger,na.rm = na_rm_setting)
  mean_rand_RandomImp <- mean(rand_RandomImp,na.rm = na_rm_setting)
  mean_rand_mice_pmm <- mean(rand_mice_pmm,na.rm = na_rm_setting)
  mean_rand_mice_cart <- mean(rand_mice_cart,na.rm = na_rm_setting)
  mean_rand_amelia <- mean(rand_amelia,na.rm = na_rm_setting)
  
  sd_rand_ClustImpute <- sd(rand_ClustImpute,na.rm = na_rm_setting)
  sd_rand_ClustImpute_false <- sd(rand_ClustImpute_false,na.rm = na_rm_setting)
  sd_rand_ClustImpute_no_weight <- sd(rand_ClustImpute_no_weight,na.rm = na_rm_setting)
  sd_rand_missRanger <- sd(rand_missRanger,na.rm = na_rm_setting)
  sd_rand_RandomImp <- sd(rand_RandomImp,na.rm = na_rm_setting)
  sd_rand_mice_pmm <- sd(rand_mice_pmm,na.rm = na_rm_setting)
  sd_rand_mice_cart <- sd(rand_mice_cart,na.rm = na_rm_setting)
  sd_rand_amelia <- sd(rand_amelia,na.rm = na_rm_setting)
  
  results_mean_randIndex[[count]] <- c(ClustImpute=mean_rand_ClustImpute,
                                  ClustImpute_assign_with_wf_false = mean_rand_ClustImpute_false,
                                  ClustImpute_no_weight=mean_rand_ClustImpute_no_weight,
                                  RandomImp=mean_rand_RandomImp,
                                  missRanger=mean_rand_missRanger,
                                  MICE_pmm=mean_rand_mice_pmm,
                                  MICE_cart=mean_rand_mice_cart,
                                  AMELIA=mean_rand_amelia)
  
  results_sd_randIndex[[count]] <- c(ClustImpute=sd_rand_ClustImpute,
                                     ClustImpute_assign_with_wf_false = sd_rand_ClustImpute_false,
                                     ClustImpute_no_weight=sd_rand_ClustImpute_no_weight,
                                     RandomImp=sd_rand_RandomImp,
                                     missRanger=sd_rand_missRanger,
                                     MICE_pmm=sd_rand_mice_pmm,
                                     MICE_cart=sd_rand_mice_cart,
                                     AMELIA=sd_rand_amelia) 
}

### Rand index
mean_randtbl <- data.frame(matrix(unlist(results_mean_randIndex),nrow=length(results_mean_randIndex), byrow=T))
colnames(mean_randtbl) <- names(results_mean_randIndex[[1]])
mean_randtbl$p <- P
mean_randtbl
write.csv(mean_randtbl,paste(type_missing,"-","Mean RandIndices Iris.csv"))

sd_randtbl <- data.frame(matrix(unlist(results_sd_randIndex),nrow=length(results_sd_randIndex), byrow=T))
colnames(sd_randtbl) <- names(results_sd_randIndex[[1]])
sd_randtbl$p <- P
sd_randtbl
write.csv(sd_randtbl,paste(type_missing,"-","SD RandIndices Iris.csv"))

### plot for rand index
data2plot <- mean_randtbl %>% pivot_longer(-p,names_to = "Algorithm",values_to = "Rand_Index")
data2plot_sd <- sd_randtbl %>% pivot_longer(-p,names_to = "Algorithm",values_to = "Rand_Index_SD")
data2plot$Rand_Index_SD <- data2plot_sd$Rand_Index_SD

pd <- position_dodge(width = 0.02)
rand_plot <- ggplot(data2plot,aes(x=p,y=Rand_Index,color=Algorithm)) + geom_line(size=1) + # + geom_point(size=2)
  theme_cowplot() +  theme(legend.position="bottom") + 
  xlab("Share of missings") + ylab("Rand Index") + scale_x_continuous(breaks = 1:7/10) + scale_y_continuous(breaks = 0:10/10) + 
  geom_pointrange(aes(ymin = Rand_Index-2*Rand_Index_SD, ymax = Rand_Index+2*Rand_Index_SD), position = pd, linetype="dotted")
rand_plot

# save results
save_plot(paste(type_missing,"-","Rand index iris.png"),rand_plot, base_width = 9, base_height = 5)

# pdf(paste(type_missing,"-","Rand index iris.pdf"),width=10)
# print(rand_plot)
# dev.off()

save(data2plot,rand_plot,file=paste(type_missing,"-","results_benchmarking_iris.Rdata"))

