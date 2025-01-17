library(ClustImpute)
library(ClusterR)
library(missForest)
library(mice)
library(Amelia)
library(missRanger)
library(microbenchmark)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(MASS)


create_random_data <- function(n,nr_other_vars=4,seedvar=739,rho=.3) {
  n <- round(n/3)*3
  set.seed(seedvar)
  
  # create correlated other variables
  sigma <- matrix(rho, nr_other_vars, nr_other_vars) + (1-rho)*diag(nr_other_vars)
  mat <- mvrnorm(n = n, mu = rep(0,nr_other_vars), Sigma = sigma)
  
  me<-4 # mean
  x <- c(rnorm(n/3,me/2,1),rnorm(2*n/3,-me/2,1)) 
  y <- c(rnorm(n/3,0,1),rnorm(n/3,me,1),rnorm(n/3,-me,1))
  true_clust <- c(rep(1,n/3),rep(2,n/3),rep(3,n/3)) # true clusters
  dat <- cbind(mat,x,y)
  dat<- as.data.frame(scale(dat)) # scaling
  return(list(data=dat,true_clusters=true_clust))
}

### Parameters
# ClustImpute
nr_iter <- 14 # iterations of procedure
n_end <- 10 # step until convergence of weight function to 1
nr_cluster <- 3 # number of clusters
c_steps <- 50 # number of cluster steps per iteration
# Clustering based on imputation
nr_iter_other <- (nr_iter-n_end) * c_steps # comparable number of steps "after" imputation

param_times <- 10 # how often to compute the benchmark

### for various n and p
N <- 2^(0:4)*100
# N <- c(50,100) # for testing
p <- .2 # .5

results <- list()
count <- 0

for (n in N) {
  count <- count + 1
  # Create random data with missings
  random_data <- create_random_data(n)
  dat <- random_data$data
  dat_with_miss <- miss_sim(dat,p,seed_nr=123)
  
  ### Benchmark time and calculate rand index
  rand_ClustImpute <- 0
  # rand_missForest <- 0
  rand_missRanger <- 0
  rand_RandomImp <- 0
  rand_mice_pmm <- 0
  rand_mice_cart <- 0
  rand_amelia <- 0
  
  set.seed(987) # for reproducibility
  mbm <- microbenchmark("ClustImpute" = {
    res <- ClustImpute(dat_with_miss,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end, seed_nr = random_seed)
    rand_ClustImpute <- rand_ClustImpute + external_validation(random_data$true_clusters, res$clusters)
  }, "RandomImp"={
    dat_random_imp <- dat_with_miss
    for (j in 1:dim(dat)[2]) {
      dat_random_imp[,j] <- impute(dat_random_imp[,j],fun="random")
    }
    cl_RandomImp <- KMeans_arma(data=dat_random_imp,clusters=3,n_iter=nr_iter_other,seed=random_seed)
    pred_RandomImp <- predict_KMeans(dat_random_imp,cl_RandomImp)
    class(pred_RandomImp) <- "numeric"
    rand_RandomImp <- rand_RandomImp + external_validation(random_data$true_clusters, pred_RandomImp)
  },
  # "MissForest"={
  #   if (n>2700) {
  #     print("MissForrest to slow")
  #     rand_missForest <- 0
  #   } else {
  #     imp_missForest <- missForest(dat_with_miss)
  #     cl_missForest <- KMeans_arma(data=imp_missForest$ximp,clusters=3,n_iter=nr_iter_other,seed=random_seed)
  #     pred_missForest <- predict_KMeans(imp_missForest$ximp,cl_missForest)
  #     class(pred_missForest) <- "numeric"
  #     rand_missForest <- rand_missForest + external_validation(random_data$true_clusters, pred_missForest)
  #   }
  # },
  "MissRanger"={
    imp_missRanger <- missRanger(dat_with_miss,pmm.k = 3)
    cl_missRanger <- KMeans_arma(data=imp_missRanger,clusters=3,n_iter=nr_iter_other,seed=random_seed)
    pred_missRanger <- predict_KMeans(imp_missRanger,cl_missRanger)
    class(pred_missRanger) <- "numeric"
    rand_missRanger <- rand_missRanger + external_validation(random_data$true_clusters, pred_missRanger)
  }, 
  "MICE_pmm"={
    dat_mice <- mice(dat_with_miss,m=1,maxit=50,meth='pmm') # single data set
    dat_mice <- complete(dat_mice)
    cl_mice <- KMeans_arma(data=dat_mice,clusters=3,n_iter=nr_iter_other,seed=random_seed)
    pred_mice <- predict_KMeans(dat_mice,cl_mice)
    class(pred_mice) <- "numeric"
    rand_mice_pmm <- rand_mice_pmm + external_validation(random_data$true_clusters, pred_mice)
  }, 
  "MICE_cart"={
    dat_mice <- mice(dat_with_miss,m=1,maxit=50,meth='cart') # single data set
    dat_mice <- complete(dat_mice)
    cl_mice <- KMeans_arma(data=dat_mice,clusters=3,n_iter=nr_iter_other,seed=random_seed)
    pred_mice <- predict_KMeans(dat_mice,cl_mice)
    class(pred_mice) <- "numeric"
    rand_mice_cart <- rand_mice_cart + external_validation(random_data$true_clusters, pred_mice)
  }, "AMELIA"= {
    # single data set, no screen output, empirical (or ridge) prior
    dat_amelia <- amelia(dat_with_miss,m=1,p2s=0,empri = .01 * nrow(dat_with_miss),autopri = 0.05) 
    dat_amelia <- dat_amelia$imputations$imp1
    cl_amelia <- KMeans_arma(data=dat_amelia,clusters=3,n_iter=nr_iter_other,seed=random_seed)
    pred_amelia <- predict_KMeans(dat_amelia,cl_amelia)
    class(pred_amelia) <- "numeric"
    rand_amelia <- rand_amelia + external_validation(random_data$true_clusters, pred_amelia)
    
  }  ,
  times = param_times, setup = {random_seed=round(runif(1)*1e5)}, unit = "s")
  
  # compute average rand index
  rand_ClustImpute <- rand_ClustImpute/param_times
  # rand_missForest <- rand_missForest/param_times
  rand_missRanger <- rand_missRanger/param_times
  rand_RandomImp <- rand_RandomImp/param_times
  rand_mice_pmm <- rand_mice_pmm/param_times
  rand_mice_cart <- rand_mice_cart/param_times
  rand_amelia <- rand_amelia/param_times
  
  results$randIndex[[count]] <- c(ClustImpute=rand_ClustImpute,RandomImp=rand_RandomImp,#missForest=rand_missForest,
                                  missRanger=rand_missRanger,MICE_pmm=rand_mice_pmm,MICE_cart=rand_mice_cart,AMELIA=rand_amelia)

  results$benchmark[[count]] <- mbm
  
  mbm_median <- print(mbm)$median
  results$benchmark_median[[count]] <- mbm_median
}

# save results
save(results,file="results_benchmarking.Rdata")


### Rand index
results$randIndex

randtbl <- data.frame(matrix(unlist(results$randIndex),nrow=length(results$randIndex), byrow=T))
colnames(randtbl) <- names(results$randIndex[[1]])
randtbl
write.csv(randtbl,"Rand indices.csv")


### Scalability plot
data2plot <- data.frame(median=unlist(results$benchmark_median))
data2plot$method <- rep(Hmisc::Cs(ClustImpute,RandomImp,MissRanger,MICE_pmm,MICE_cart,AMELIA),length(N))
data2plot$n <- rep(N, each=length(results$randIndex[[1]]))

# with shared legend
ps1 <- ggplot(data2plot,aes(x=n,y=median,color=method)) + geom_point() + theme_cowplot() + geom_smooth() + 
  xlab("Numer of observations n") + ylab("Median running time in seconds")
legend <- get_legend(ps1)
ps2 <- ggplot(data2plot,aes(x=n,y=median,color=method)) + geom_point() + theme_cowplot() + geom_smooth() + 
  scale_y_log10() + xlab("Numer of observations n") + ylab("Median running time in seconds") + theme(legend.position="none")
p2 <- plot_grid(ps1 + theme(legend.position="none"),ps2,legend,nrow=1,rel_widths = c(3,3,1))
p2
save_plot("Scalability.png",p2,base_aspect_ratio = 2.5)
