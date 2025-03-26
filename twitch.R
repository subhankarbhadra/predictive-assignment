source("predinf_dcbm.R")

#Clustering algorithm (SVD + K-means) (regularized)
fast.clustering_DCBM_reg <- function(A, k, niter, nstart, reg = "none") {
  if(reg == "none") {
    t1 <- system.time(e <- irlba(A, nu = k, nv = k))[3]
  } else if(reg == "laplacian") {
    M <- -laplacian_matrix(graph_from_adjacency_matrix(A, mode = "undirected"), normalized = T)
    diag(M) <- 0
    
    t1 <- system.time(e <- irlba(M, nu = k, nv = k))[3]
  } else if(reg == "laplacian+deg") {
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

#Spectral clustering on full matrix
specclustering_DCBM_reg <- function(A, k, niter, nstart, reg) {
  n <- ncol(A)
  deg <- colSums(A)
  
  #excluding 0-degree nodes from the graph
  d0 <- which(deg != 0)
  d1 <- setdiff(1:n, d0)
  cw <- fast.clustering_DCBM_reg(A[d0, d0], k, niter, nstart, reg)
  pred <- integer(n)
  pred[d0] <- cw$cluster
  pred[d1] <- sample(1:k, length(d1), T)
  list(cluster = pred, times = cw$times) 
}

#Clustering based on predictive assignment
effclustering_DCBM_reg <- function(A, k, p, rw = F, method = "np", niter = niter, nstart = nstart, reg) {
  #sampling
  t1 <- system.time({
    n <- ncol(A);
    m <- floor(n^p);
    if(rw) {
      subnet <- crwsample(A@i, A@p, m)
    } else {
      subnet <- sample(1:n, m)
    };
    
    #excluding 0-degree nodes from the subgraph
    zd <- which(colSums(A[subnet, subnet]) == 0);
    if(length(zd) > 0) subnet <- subnet[-zd];
    rest <- setdiff(1:n, subnet);
    m <- length(subnet)
  })[3]
  
  #spectral clustering on subgraph
  t2 <- system.time(cw <- fast.clustering_DCBM_reg(A[subnet, subnet], k, niter, nstart, reg))[3]
  
  #predictive assignment
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
eval_eff <- function(p, reg){
  time <- system.time(out <- effclustering_DCBM_reg(A, k, p, T, "np", niter, nstart, reg))[3]
  err <- calc_error_eff(out, classes)
  c(err, out$times, unname(time))
}
dd <- rbind(cbind(t(replicate(N, eval_eff(0.8, "none"))), rep(0.8, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.85, "none"))), rep(0.85, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.9, "none"))), rep(0.9, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.8, "laplacian+deg"))), rep(0.8, N), rep("laplacian+deg", N)),
            cbind(t(replicate(N, eval_eff(0.85, "laplacian+deg"))), rep(0.85, N), rep("laplacian+deg", N)),
            cbind(t(replicate(N, eval_eff(0.9, "laplacian+deg"))), rep(0.9, N), rep("laplacian+deg", N)))

dd <- as.data.frame(dd[,c(1:3,7:9)])
dd$V1 <- as.numeric(dd$V1)
dd$V2 <- as.numeric(dd$V2)
dd$V3 <- as.numeric(dd$V3)
dd$V4 <- as.numeric(dd$V4)
tb <- dd %>% group_by(V6, V5) %>% 
  summarise(avg.rate1 = mean(V1)*100,
            sd.rate1 = sd(V1)*100,avg.rate2 = mean(V2)*100,
            sd.rate2 = sd(V2)*100,avg.rate = mean(V3)*100,
            sd.rate = sd(V3)*100, avg.tot = mean(V4)) %>% as_tibble()
tb
tb$avg.rate

eval_spec <- function(){
  time <- system.time(out <- specclustering_DCBM_reg(A, k, niter, nstart, "none"))[3]
  err <- calc_error_spec(out, classes)
  c(err, out$times, unname(time))
}
eval_spec()[c(1,5)]

eval_spec <- function(){
  time <- system.time(out <- specclustering_DCBM_reg(A, k, niter, nstart, "laplacian+deg"))[3]
  err <- calc_error_spec(out, classes)
  c(err, out$times, unname(time))
}

eval_spec()[c(1,5)]
