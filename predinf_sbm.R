library(tidyverse)
library(parallel)
library(igraph)
library(Matrix)
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

# Performs spectral clustering
#
# Arguments:
#  M: n x n symmetric matrix
#  k: Integer, number of communities
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd and kmeans
#
fast.clustering <- function(M, k, niter, nstart) {
  t1 <- system.time(e <- irlba(M, nu = k, nv = k))[3]
  S <- e$u
  t2 <- system.time(c <- kmeans(S, k, iter.max = niter, nstart = nstart)$cluster)[3]
  list(cluster = c, times = c(t1, t2))
}

# Performs spectral clustering given an adjacency matrix
# Same as `fast.clustering`, with an indicator for whether it is applied on the adjacency matrix itself, or its bias-adjusted version 
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  adjusted: Indicator for whether spectral clustering is to be applied on A, or its bias-adjusted version, A^2 - D
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd and kmeans
#
specclustering <- function(A, k, adjusted = F, niter, nstart) {
  n <- ncol(A)
  if(adjusted == F) M <- A
  else {
    M <- as.matrix(crossprod(A))
    M <- `diag<-`(M, 0)
  }
  cw <- fast.clustering(M, k, niter, nstart)
  cw
}


# Performs Step 3 of predictive assignment using closest community approach
# This is a computationally efficient implementation using `inline`
# Arguments:
#  x: (n-m) x k numeric matrix, the estimated theta matrix multiplied by (-2)
#  l: list of length (n-m), neighbors of the remaining nodes
#  v: numeric vector of length k,  squared l_2 norms of columns of the estimated theta matrix 
# Returns:
#  Integer vector, estimated community assignments
#
sig <- c(x = "double", l = "list", v = "double")

bod <- '
  double *px = REAL(x);
  double *pv = REAL(v);
  R_xlen_t nx = XLENGTH(x);
  int *d = INTEGER(getAttrib(x, R_DimSymbol));
  int m = d[0];
  int n = d[1];
  R_xlen_t N = XLENGTH(l);
  SEXP res = PROTECT(allocVector(INTSXP, N));
  int *pres = INTEGER(res);
  SEXP index;
  R_xlen_t nindex;
  int *pindex;
  double sum;
  int min_index;
  double min_sum;

  for (R_xlen_t i = 0; i < N; ++i)
  {
    index = VECTOR_ELT(l, i);
    nindex = XLENGTH(index);
    pindex = INTEGER(index);
    min_index = -1;
    min_sum = m;
  for (R_xlen_t xpos = 0, colpos = 0; xpos < nx; xpos += m, ++colpos)
  {
    sum = 0.0;
    for (R_xlen_t k = 0; k < nindex; ++k)
    {
      sum += px[xpos + pindex[k] - 1];
    }
    sum += pv[colpos];

    // Update min_index and min_sum if the current sum is smaller
    if (sum < min_sum) {
      min_sum = sum;
      min_index = colpos; // Update index with the minimum sum
    }
  }
  pres[i] = min_index + 1;
}
  UNPROTECT(1);
  return res;
  '

min_dist_inline <- cfunction(sig, bod, language = "C")


# Clustering based on predictive assignment (closest community approach)
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  p: Numeric, (log m)/(log n) where m is the size of the subgraph
#  method: Predictive assignment approach ("cc" for closest community)
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  subsample: Integer vector, nodes selected in the subgraph
#  times: runtimes for sampling, spectral clustering on subgraph and predictive assignment of remaining nodes
#
effclustering <- function(A, k, p, method = "cc", adjusted = F, niter, nstart) {
  # sampling
  t1 <- system.time({
    n <- ncol(A);
    m <- floor(n^p);
    subnet <- sample(1:n, m);
    rest <- setdiff(1:n, subnet);
    if(adjusted == F) M <- A[subnet, subnet]
    else {
      M <- as.matrix(tcrossprod(A[subnet,]))
      M <- `diag<-`(M, 0)
    }
  })[3]
  
  # community detection on subgraph
  t2 <- system.time(cw <- fast.clustering(M, k, niter, nstart))[3]
  
  # predictive assignment of remaining nodes
  t3 <- system.time({
    mr <- unname(table(cw$cluster));
    if(method == "cc") {
      # estimate theta
      theta <- sapply(1:k, function(r) rowSums(A[rest, subnet[cw$cluster == r], drop = F])/mr[r])
      
      # find neighbors of nodes
      Arest <- A[rest, rest]
      a <- Arest@i + 1L;
      b <- Arest@p + 1L;
      nbr <- list();
      for (i in seq_along(b[-1])) {
        nbr[[i]] <- a[b[i]:(b[i+1]-1)]
      }
      l <- length(nbr)
      if(l != (n-m)) stop("neighbor extraction error")
      
      # find community with minimum distance
      tk <- colSums(theta^2)
      restc <- min_dist_inline(-2*theta, nbr, tk)
    } else stop("wrong method")
  })[3]
  
  pred <- numeric(n)
  pred[subnet] <- cw$cluster
  pred[rest] <- restc
  
  list(cluster = pred, subsample = subnet, times = c(t1, t2, t3))
}

# Calculates misclassification error from spectral clustering results
# Arguments:
#  out: An output from specclustering()
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
#  out: An output from effclustering()
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