source("predinf_dcbm.R")

# Performs regularized spectral clustering with or without laplacian transformation
#
# Arguments:
#  A: n x n symmetric binary adjacency matrix
#  k: Integer, number of communities
#  niter: Integer, maximum number of iterations for kmeans
#  nstart: Integer, number of initializations for kmeans
#  version: "none" for without laplacian, "laplacian" for with laplacian, "laplacian+deg" for with laplacian and degree-regularization
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd, row-normalization and kmeans
#
fast.clustering_DCBM_reg <- function(A, k, niter, nstart, version = "none") {
  if(version == "none") {
    t1 <- system.time(e <- irlba(A, nu = k, nv = k))[3]
  } else if(version == "laplacian") {
    M <- -laplacian_matrix(graph_from_adjacency_matrix(A, mode = "undirected"), normalized = T)
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
#  version: "none" for without laplacian, "laplacian" for with laplacian, "laplacian+deg" for with laplacian and degree-regularization
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  times: runtimes for performing svd and kmeans
#
specclustering_DCBM_reg <- function(A, k, niter, nstart, version) {
  n <- ncol(A)
  deg <- colSums(A)
  
  #excluding 0-degree nodes from the graph
  d0 <- which(deg != 0)
  d1 <- setdiff(1:n, d0)
  cw <- fast.clustering_DCBM_reg(A[d0, d0], k, niter, nstart, version)
  pred <- integer(n)
  pred[d0] <- cw$cluster
  pred[d1] <- sample(1:k, length(d1), T)
  list(cluster = pred, times = cw$times) 
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
#  version: "none" for without laplacian, "laplacian" for with laplacian, "laplacian+deg" for with laplacian and degree-regularization
#
# Returns: A list
#  cluster: Integer vector, estimated community assignments
#  subsample: Integer vector, nodes selected in the subgraph
#  times: runtimes for sampling, spectral clustering on subgraph and predictive assignment of remaining nodes
#
effclustering_DCBM_reg <- function(A, k, p, rw = F, method = "np", niter = niter, nstart = nstart, version) {
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
  t2 <- system.time(cw <- fast.clustering_DCBM_reg(A[subnet, subnet], k, niter, nstart, version))[3]
  
  # predictive assignment of remaining nodes
  t3 <- system.time(if(method == "np") {
    #excluding nodes with zero connection to the subgraph 
    d0 <- which(colSums(A[subnet, rest]) == 0)
    rest0 <- rest[d0]
    rest1 <- setdiff(rest, rest0)
    deg.rest1 <- colSums(A[subnet, rest1])
    
    comm_sub <- split(subnet, cw$cluster)
    
    #calculating a scaled version of omega-hat
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

sourceCpp("crwsample.cpp")

N <- 100; niter <- 1000; nstart <- 2000
twitch <- readRDS("twitch.RData")
A <- twitch$adjacency
classes <- twitch$classes
k <- length(unique(classes))

set.seed(12345)

# Applies predictive assignment on twitch network for different values of p = (log m)/(log n)
# Returns community detection errors for subgraph nodes, remaining nodes, and all nodes, and average runtime
#
eval_eff <- function(p, version){
  time <- system.time(out <- effclustering_DCBM_reg(A, k, p, T, "np", niter, nstart, version))[3]
  err <- calc_error_eff(out, classes)
  c(err, unname(time))
}
dd <- rbind(cbind(t(replicate(N, eval_eff(0.8, "none"))), rep(0.8, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.85, "none"))), rep(0.85, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.9, "none"))), rep(0.9, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.8, "laplacian+deg"))), rep(0.8, N), rep("laplacian+deg", N)),
            cbind(t(replicate(N, eval_eff(0.85, "laplacian+deg"))), rep(0.85, N), rep("laplacian+deg", N)),
            cbind(t(replicate(N, eval_eff(0.9, "laplacian+deg"))), rep(0.9, N), rep("laplacian+deg", N)))

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

# Applies spectral clustering on twitch network
eval_spec <- function(version){
  time <- system.time(out <- specclustering_DCBM_reg(A, k, niter, nstart, version))[3]
  err <- calc_error_spec(out, classes)
  c(err, unname(time))
}

# Community detection error and average runtime
eval_spec("none")
eval_spec("laplacian+deg")
