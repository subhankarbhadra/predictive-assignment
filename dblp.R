source("predinf_sbm.R")

# Performs spectral clustering with or without laplacian transformation
#
# Arguments:
#  M: n x n symmetric matrix
#  k: Integer, number of communities
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#  version: "none" for without laplacian, "laplacian" for with laplacian, "laplacian+deg" for with laplacian and degree-regularization
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd and kmeans
#
fast.clustering_reg <- function(A, k, niter, nstart, version = "none") {
  if(version == "none") {
    t1 <- system.time(e <- irlba(A, nu = k, nv = k))[3]
  } else if(version == "laplacian") {
    M <- -laplacian_matrix(graph_from_adjacency_matrix(A, mode = "undirected"), normalization = "symmetric")
    diag(M) <- 0
    
    t1 <- system.time(e <- irlba(M, nu = k, nv = k))[3]
  } else if(version == "laplacian+deg") {
    n <- ncol(A)
    degree <- colSums(A)
    
    D <- sparseMatrix(i = 1:n, j = 1:n, x = 1/sqrt(degree + mean(degree)))
    M <- tcrossprod(crossprod(D, A), D)
    
    t1 <- system.time(e <- irlba(M, nu = k, nv = k))[3]
  }
  
  S <- e$u
  t2 <- system.time(c <- kmeans(S, k, iter.max = niter, nstart = nstart)$cluster)[3]
  list(cluster = c, times = c(t1, t2))
}

# Performs spectral clustering given an adjacency matrix with or without laplacian transformation
# Same as `fast.clustering_reg`, with an indicator for whether it is applied on the adjacency matrix itself, or its bias-adjusted version 
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  adjusted: Indicator for whether spectral clustering is to be applied on A, or its bias-adjusted version, A^2 - D
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#  version: "none" for without laplacian, "laplacian" for with laplacian, "laplacian+deg" for with laplacian and degree-regularization
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd and kmeans
#
specclustering_reg <- function(A, k, adjusted = F, niter, nstart, version) {
  n <- ncol(A)
  if(adjusted == F) M <- A
  else {
    M <- as.matrix(crossprod(A))
    M <- `diag<-`(M, 0)
  }
  cw <- fast.clustering_reg(M, k, niter, nstart, version)
  cw
}

# Clustering based on predictive assignment (closest community approach)
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  p: Numeric, (log m)/(log n) where m is the size of the subgraph
#  method: Predictive assignment approach ("cc" for closest community)
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#  version: "none" for without laplacian, "laplacian" for with laplacian, "laplacian+deg" for with laplacian and degree-regularization
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  subsample: Integer vector, nodes selected in the subgraph
#  times: runtimes for sampling, spectral clustering on subgraph and predictive assignment of remaining nodes
#
effclustering_reg <- function(A, k, p, method = "cc", adjusted = F, niter, nstart, version = "none") {
  #sampling
  t1 <- system.time({
    n <- ncol(A);
    m <- ceiling(n^p);
    subnet <- sample(1:n, m);
    rest <- setdiff(1:n, subnet);
    if(adjusted == F) M <- A[subnet, subnet]
    else {
      M <- as.matrix(tcrossprod(A[subnet,]))
      M <- `diag<-`(M, 0)
    }
  })[3]
  
  #spectral clustering on subgraph
  t2 <- system.time(cw <- fast.clustering_reg(M, k, niter, nstart, version))[3]
  
  #predictive assignment of remaining nodes
  t3 <- system.time({
    mr <- unname(table(cw$cluster));
    if(method == "cc") {
      #calculate theta
      theta <- sapply(1:k, function(r) rowSums(A[rest, subnet[cw$cluster == r], drop = F])/mr[r])
      
      #find neighbors of nodes
      Arest <- A[rest, rest]
      a <- Arest@i + 1L;
      b <- Arest@p + 1L;
      nbr <- list();
      for (i in seq_along(b[-1])) {
        nbr[[i]] <- a[b[i]:(b[i+1]-1)]
      }
      l <- length(nbr)
      if(l != (n-m)) stop("neighbor extraction error")
      
      #find community with minimum distance
      tk <- colSums(theta^2)
      restc <- min_dist_inline(-2*theta, nbr, tk)
    } else stop("wrong method")
  })[3]
  
  pred <- numeric(n)
  pred[subnet] <- cw$cluster
  pred[rest] <- restc
  
  list(cluster = pred, subsample = subnet, times = c(t1, t2, t3))
}

N <- 100; niter <- 1000; nstart <- 2000
dblp <- readRDS("dblp.RData")
A <- dblp$adjacency
classes <- dblp$classes
k <- length(unique(classes))

set.seed(12345)

# Applies predictive assignment on DBLP network for different values of p = (log m)/(log n)
# Returns community detection errors for subgraph nodes, remaining nodes, and all nodes, and average runtime
#
eval_eff <- function(p, version){
  time <- system.time(out <- effclustering_reg(A, k, p, "cc", adjusted = F, niter, nstart, version))[3]
  err <- calc_error_eff(out, classes)
  c(err, unname(time))
}

dd <- rbind(cbind(t(replicate(N, eval_eff(0.70, "none"))), rep(0.70, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.75, "none"))), rep(0.75, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.80, "none"))), rep(0.80, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.70, "laplacian"))), rep(0.70, N), rep("laplacian", N)),
            cbind(t(replicate(N, eval_eff(0.75, "laplacian"))), rep(0.75, N), rep("laplacian", N)),
            cbind(t(replicate(N, eval_eff(0.80, "laplacian"))), rep(0.80, N), rep("laplacian", N)))

dd <- as.data.frame(dd)
colnames(dd) <- c("rate1", "rate2", "rate", "time", "p", "version")

# Summarizing results 
# Average community detection error and average runtime
#
for (i in 1:4) dd[[i]] <- as.numeric(dd[[i]])
tb <- dd %>% group_by(version, p) %>% 
  summarise(avg.rate = mean(rate)*100, sd.rate = sd(rate)*100, 
            avg.time = mean(time)) %>% as_tibble()
tb

# Applies spectral clustering on DBLP network
eval_spec <- function(version){
  time <- system.time(out <- specclustering_reg(A, k, F, niter, nstart, version))[3]
  err <- calc_error_spec(out, classes)
  c(err, unname(time))
}

# Community detection error and average runtime
eval_spec("none")
eval_spec("laplacian")
