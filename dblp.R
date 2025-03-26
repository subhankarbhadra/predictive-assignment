source("predinf_sbm.R")

#Clustering algorithm (SVD + K-means) (regularized)
fast.clustering_reg <- function(A, k, niter, nstart, reg = "none") {
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
  t2 <- system.time(c <- kmeans(S, k, iter.max = niter, nstart = nstart)$cluster)[3]
  list(cluster = c, times = c(t1, t2))
}

#Spectral clustering on full matrix (regularized)
specclustering_reg <- function(A, k, adjusted = F, niter, nstart, reg) {
  n <- ncol(A)
  if(adjusted == F) M <- A
  else {
    M <- as.matrix(crossprod(A))
    M <- `diag<-`(M, 0)
  }
  cw <- fast.clustering_reg(M, k, niter, nstart, reg)
  cw
}

#Clustering based on predictive assignment
effclustering_reg <- function(A, k, p, method = "cc", adjusted = F, niter, nstart, reg = "none") {
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
  t2 <- system.time(cw <- fast.clustering_reg(M, k, niter, nstart, reg))[3]
  
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

eval_eff <- function(p, reg){
  time <- system.time(out <- effclustering_reg(A, k, p, "cc", adjusted = F, niter, nstart, reg))[3]
  err <- calc_error_eff(out, classes)
  c(err, out$times, unname(time))
}

dd <- rbind(cbind(t(replicate(N, eval_eff(0.70, "none"))), rep(0.70, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.75, "none"))), rep(0.75, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.80, "none"))), rep(0.80, N), rep("none", N)),
            cbind(t(replicate(N, eval_eff(0.70, "laplacian"))), rep(0.70, N), rep("laplacian", N)),
            cbind(t(replicate(N, eval_eff(0.75, "laplacian"))), rep(0.75, N), rep("laplacian", N)),
            cbind(t(replicate(N, eval_eff(0.80, "laplacian"))), rep(0.80, N), rep("laplacian", N)))

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

eval_spec <- function(reg){
  time <- system.time(out <- specclustering_reg(A, k, F, niter, nstart, reg))[3]
  err <- calc_error_spec(out, classes)
  c(err, out$times, unname(time))
}
eval_spec("none")[c(1,4)]; eval_spec("laplacian")[c(1,4)]
