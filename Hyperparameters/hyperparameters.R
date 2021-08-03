#### Init ####
library(tidyverse)
library(ClustImpute)
library(ClusterR)
library(missRanger)
library(microbenchmark)
library(Hmisc)
library(cowplot)
library(RColorBrewer)

#### Data prep ####

data("iris")
dat <- as.data.frame(scale(iris[,1:4]))
true_label <- as.numeric(iris$Species)

type_missing <- "MAR" # use MCAR or MAR

# how often to run each configuration
nr_runs <- 30

### Parameters
nr_iter <- c(5,10,20,30,50) # iterations of procedure
nr_cluster <- 3 # number of clusters
c_steps <- c(1,10,25,50) # numer of cluster steps per iteration
# make grid
ClustImpute_params <- expand_grid(nr_iter,c_steps)
# no weight function
ClustImpute_params$n_end <- 1 # no weight function used
ClustImpute_params$group <- "no weight function"
# weight function, varying convergence
n_end_frac <- c(.3,.6,.9) # define n_end as a fraction of nr_iter
ClustImpute_params_base <-  ClustImpute_params
for (j in 1:length(n_end_frac)) {
  tmp <- ClustImpute_params_base
  tmp$n_end <- round(n_end_frac[j]*ClustImpute_params_base$nr_iter)
  tmp$group <- paste("WF converges at ",n_end_frac[j])
  ClustImpute_params <- rbind(ClustImpute_params,tmp)
}

ClustImpute_params$all_iter <- ClustImpute_params$nr_iter * ClustImpute_params$c_steps

### share of missings
p <- .3

set.seed(124)
# create data with missings
dat_with_miss <- miss_sim(dat,p,seed_nr=739,type=type_missing)


#### Runs ####

for (k in 1:nrow(ClustImpute_params)) {
  rand_ClustImpute <- 0
  for (runs in 1:nr_runs) {
    res <- ClustImpute(dat_with_miss, nr_cluster = nr_cluster, 
                       nr_iter=as.integer(ClustImpute_params[1,"nr_iter"]), 
                       c_steps=as.integer(ClustImpute_params[k,"c_steps"]), 
                       n_end=as.integer(ClustImpute_params[k,"n_end"]),
                       seed_nr = 125+k+runs)
    rand_ClustImpute <- rand_ClustImpute + external_validation(true_label, res$clusters)
  }
  ClustImpute_params[k,"value"] <- rand_ClustImpute/nr_runs
}

#### Save ####

save(ClustImpute_params,file="ClustImpute_params.Rdata")

#### Viz ####

plt <- ggplot(ClustImpute_params,aes(x=nr_iter,y=value,color=factor(c_steps),size=all_iter)) + geom_point() + facet_grid(~group) + scale_color_brewer(palette = "Dark2")

# plotly::ggplotly(plt)

plt_journal <- plt + theme_cowplot() + xlab("Number of iterations: nr_iter") + ylab("Rand Index") + scale_y_continuous(breaks = seq(46,60,2)/100) +
  scale_x_continuous(breaks = unique(ClustImpute_params$nr_iter)) +  labs(color='c_steps') + labs(size="nr_iter * c_steps")
plt_journal

save_plot("Result hyper parameters.png",plt_journal,base_aspect_ratio = 2.5)

#### Descriptive ####

ClustImpute_params %>% arrange(-value)

ClustImpute_params %>% group_by(nr_iter) %>% summarise(mean(value))
ClustImpute_params %>% group_by(group) %>% summarise(mean(value))
ClustImpute_params %>% group_by(c_steps) %>% summarise(mean(value))

ClustImpute_params %>% group_by(nr_iter) %>% summarise(max(value))
ClustImpute_params %>% group_by(group) %>% summarise(max(value))
ClustImpute_params %>% group_by(c_steps) %>% summarise(max(value))
