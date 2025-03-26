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

#Calculate errors from predictive assignment using SC for different values of p = (log m)/(log n)
ratecal <- function(networks, h, k, classes, niter, nstart) {
  N <- length(networks)
  
  rate1 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.85, "cc", F, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate2 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.90, "cc", F, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate3 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering(networks[[i]], k, 0.95, "cc", F, niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  
  
  d1 <- as.data.frame(cbind(rbind(rate1, rate2, rate3), rep("cc", N*3), rep("SC", N*3), rep(h, N*3), rep(c(0.85, 0.90, 0.95), each = N)))
  names(d1) <- c("rate1", "rate2", "rate", "t1", "t2", "t3", "time", "method", "comdet", "h", "p")
  
  rate1 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- specclustering(networks[[i]], k, F, niter, nstart))[3]
    err <- calc_error_spec(out, classes)
    c(err, 0, err, 0, sum(out$times), 0, unname(time), "full", "SC")
  }))
  
  d3 <- as.data.frame(cbind(rbind(rate1), rep(h,N), rep(1, N)))
  names(d3) <- c("rate1", "rate2", "rate", "t1", "t2", "t3", "time", "method", "comdet", "h", "p")
  
  rbind(d1, d3)
}

out <- rategen_check(h = 3, n = 100000, k = 20, delta = 0.01, N = 4, niter = 1000, nstart = 2000)
out
