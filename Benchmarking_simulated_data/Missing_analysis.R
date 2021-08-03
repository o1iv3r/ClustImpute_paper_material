set.seed(739)
dim <- dim(dat)[2]

type <- "MAR"

if(type == "MAR") {
  cor_param <- 2 * stats::runif(dim * (dim - 1)/2) - 1
  cor_matrix <- stats::cov2cor(copula::p2P(cor_param) %*% 
                                 t(copula::p2P(cor_param)))
  cor_matrix <- copula::P2p(cor_matrix)
} 
if(type == "MCAR") {
  cor_matrix <- diag(dim)
  cor_matrix <- copula::P2p(cor_matrix)
}

myCop <- copula::normalCopula(param = cor_matrix, dim, dispstr = "un")
myMvd <- copula::mvdc(copula = myCop, margins = c("binom"), 
                      marginsIdentical = T, paramMargins = list(list(size = 1, 
                                                                     prob = p)))
mis_ind <- copula::rMvdc(dim(dat)[1], myMvd)
# mis_ind_NA <- sapply(mis_ind, function(x) ifelse(x == 1, 
#                                                  NA, x))
# return(dat + mis_ind_NA)

mis_ind <- as.data.frame(mis_ind)

# Correlation plot
# psych::corPlot(mis_ind,numbers = TRUE)

names(mis_ind) <- names(iris)[1:4]

library(corrplot)
corrplot(cor(mis_ind),method="number")
