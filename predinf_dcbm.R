library(tidyverse)
library(igraph)
library(Matrix)
library(parallel)
library(peakRAM)
library(irlba)
library(Rcpp)
library(inline)

# Rearranges the rows and columns of a matrix so that its largest element is moved to the (1,1) position
#
# Arguments:
#   M : Numeric matrix
#
# Returns:
#   Numeric matrix with the largest element moved to the (1,1) position
#
maxswap <- function(M) {
  n <- nrow(M)
  k <- which.max(M)
  j <- ceiling(k/n)
  i <- k - n*(j-1)
  c(i, j)
  M[c(1, i),] <- M[c(i, 1),]
  M[,c(1, j)] <- M[,c(j, 1)]
  M
}

# Greedy algorithm to correct label permutation
# Iteratively permute the rows and columns of a confusion matrix to align predicted labels with true labels
#
# Arguments:
#   M : Confusion matrix (numeric)
#
# Returns:
#   Permuted confusion matrix with improved alignment of labels
#
miscorrect <- function(M) {
  n <- nrow(M)
  for (i in 1:(n-1)) M[i:n, i:n] <- maxswap(M[i:n, i:n])
  M
}

# Performs regularized spectral clustering
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd, row-normalization and kmeans
#
fast.clustering_DCBM <- function(A, k, niter, nstart) {
  t1 <- system.time(e <- irlba(A, nu = k, nv = k))[3]
  S <- e$u
  t2 <- system.time(S_reg <- t(sapply(1:nrow(S), function(i) S[i,]/sqrt(sum(S[i,]^2)))))[3]
  t3 <- system.time(c <- kmeans(S_reg, k, iter.max = niter, nstart = nstart)$cluster)[3]
  
  list(cluster = c, times = c(t1, t2, t3))
}

# Performs regularized spectral clustering after excluding zero-degrees nodes
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd and kmeans
#
specclustering_DCBM <- function(A, k, niter, nstart) {
  n <- ncol(A)
  deg <- colSums(A)
  
  #excluding zero-degree nodes from the graph
  d0 <- which(deg != 0)
  d1 <- setdiff(1:n, d0)
  cw <- fast.clustering_DCBM(A[d0, d0], k, niter, nstart)
  pred <- integer(n)
  pred[d0] <- cw$cluster
  pred[d1] <- sample(1:k, length(d1), T)
  list(cluster = pred, times = cw$times) 
}

# Generates a network from a DCBM
#
# Arguments:
#  n: Integer, number of nodes
#  k: Integer, number of communities
#  W: k x k numeric matrix, block probability matrix
#  t: Numeric vector of degree parameters
#  m: Integer vector of true communities of nodes
#
# Returns: 
#  n x n symmetric binary adjacency matrix
#
DCBM.fast <- function(n, k, W, t, comm)
{
  on.exit(gc())
  
  stor <- lapply(1:(n-1), function(i) {
    tmp <- rbinom(n-i, 1, pmin(W[comm[i],comm[(i+1):n]]*t[i]*t[(i+1):n], 1))
    i + which(tmp == 1)
  })
  
  size <- sapply(1:(n-1), function(i) length(stor[[i]]))
  vec <- list('i' = rep(1:(n-1), size), 'j' = unlist(stor))
  
  A <- sparseMatrix(i = c(vec$i, vec$j), j = c(vec$j, vec$i), x = 1, dims = c(n,n))
}

# Clustering based on predictive assignment (node popularity approach)
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  p: Numeric, (log m)/(log n) where m is the size of the subgraph
#  rw: Sampling method (FALSE for Simple Random Sampling, TRUE for Random Walk Sampling)
#  method: Predictive assignment approach ("np" for node popularity)
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  subsample: Integer vector, nodes selected in the subgraph
#  times: runtimes for sampling, spectral clustering on subgraph and predictive assignment of remaining nodes
#
effclustering_DCBM <- function(A, k, p, rw = F, method = "np", niter = niter, nstart = nstart) {
  # sampling
  t1 <- system.time({
    n <- ncol(A);
    m <- floor(n^p);
    if(rw) {
      subnet <- crwsample(A@i, A@p, m)
    } else {
      subnet <- sample(1:n, m)
    };
    
    # excluding 0-degree nodes from the subgraph
    zd <- which(colSums(A[subnet, subnet]) == 0);
    if(length(zd) > 0) subnet <- subnet[-zd];
    rest <- setdiff(1:n, subnet);
    m <- length(subnet)
  })[3]
  
  # spectral clustering on subgraph
  t2 <- system.time(cw <- fast.clustering_DCBM(A[subnet, subnet], k, niter, nstart))[3]
  
  # predictive assignment of remaining nodes
  t3 <- system.time(if(method == "np") {
    # excluding nodes with zero connection to the subgraph 
    d0 <- which(colSums(A[subnet, rest]) == 0)
    rest0 <- rest[d0]
    rest1 <- setdiff(rest, rest0)
    deg.rest1 <- colSums(A[subnet, rest1])
    
    comm_sub <- split(subnet, cw$cluster)
    
    # calculating a scaled version of omega-hat
    What <- matrix(0, ncol = k, nrow = k)
    
    mr <- unname(table(cw$cluster))
    A_rowsum <- sapply(1:k, function(r) rowSums(A[, subnet[cw$cluster == r], drop = F]))
    for(r in 1:k) {
      What[r,] <- colSums(A_rowsum[comm_sub[[r]], ])/mr
    }
    
    Wmod <- What/rep(colSums(What), each = k)
    
    restc1 <- sapply(seq_along(rest1), function(i) {
      which.min(rowSums((t(Wmod) - A_rowsum[rep(rest1[i], k), ]/deg.rest1[i])^2))
    })
    restc0 <- sample(1:k, length(rest0), replace = T)
  } else stop("wrong method"))[3]
  
  pred <- numeric(n)
  pred[subnet] <- cw$cluster
  pred[rest1] <- restc1
  pred[rest0] <- restc0
  
  list(cluster = pred, subsample = subnet, times = c(t1, t2, t3))
}

# Calculates misclassification error from spectral clustering results
# Arguments:
#  out: An output from specclustering_DCBM()
#  classes: Integer vector, true community labels
#
# Returns:
#  Numeric, misclassification error
#
calc_error_spec <- function(out, classes) {
  pred <- out$cluster
  n <- length(pred)
  
  k <- length(unique(classes))
  pred <- factor(pred, 1:k)
  classes <- factor(classes, 1:k)
  
  cf <- table(classes, pred)
  cf <- miscorrect(cf)
  
  1 - sum(diag(cf))/n
}

# Calculates misclassification error from predictive assignment-based clustering results
# Arguments:
#  out: An output from effclustering_DCBM()
#  classes: Integer vector, true community labels
#
# Returns:
#  Numeric vector, misclassification errors for subgraph nodes, remaining nodes, and all nodes 
#
calc_error_eff <- function(out, classes) {
  pred <- out$cluster
  subnet <- out$subsample
  n <- length(pred)
  m <- length(subnet)
  rest <- setdiff(1:n, subnet)
  
  k <- length(unique(classes))
  pred <- factor(pred, 1:k)
  classes <- factor(classes, 1:k)
  
  cf1 <- table(classes[subnet], pred[subnet])
  cf2 <- table(classes[rest], pred[rest])
  cf <- table(classes, pred)
  
  cf1 <- miscorrect(cf1)
  cf2 <- miscorrect(cf2)
  cf <- miscorrect(cf)
  
  c(1 - sum(diag(cf1))/m, 1 - sum(diag(cf2))/(n-m), 1 - sum(diag(cf))/n)
}