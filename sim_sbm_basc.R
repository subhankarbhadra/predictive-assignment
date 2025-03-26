source("predinf_sbm.R")

#Simulate N networks from the SBM and apply predictive assignment
#n: number of nodes
#k: number of communities
#h: homophily factor
#delta: average density
rategen_check <- function(h, n, k, delta, N, niter, nstart) {
  classes <- sort(rep(1:k, ceiling(n/k))[1:n])
  perm <- sample.int(n, n)
  classes <- classes[perm]
  
  W <- (delta*k/(h + k - 1))*(diag(rep((h-1), k)) + outer(rep(1,k), (rep(1,k))))
  
  networks <- list(0)
  for (i in 1:N) networks[[i]] <- as_adj(sample_sbm(n, W, as.numeric(table(classes))))[perm, perm]
  
  result <- ratecal(networks, h, k, classes, niter, nstart)
  result
}

#Calculate errors from predictive assignment using BASC for different values of p = (log m)/(log n)
ratecal <- function(networks, h, k, classes, niter, nstart) {
  N <- length(networks)
  
  rate1 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.70, "cc", T, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate2 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.75, "cc", T, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate3 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.80, "cc", T, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate4 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.85, "cc", T, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  
  d1 <- as.data.frame(cbind(rbind(rate1, rate2, rate3, rate4), rep("cc", N*4), rep("BASC", N*4), rep(h, N*4), rep(c(0.70,0.75,0.80,0.85), each = N)))
  names(d1) <- c("rate1", "rate2", "rate", "t1", "t2", "t3", "time", "method", "comdet", "h", "p")
  
  d1
}

out <- rategen_check(h = 3, n = 100000, k = 20, delta = 0.01, N = 4, niter = 1000, nstart = 2000)
out
