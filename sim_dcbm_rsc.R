source("predinf_dcbm.R")

#Simulate N networks from the DCBM and apply predictive assignment
#n: number of nodes
#k: number of communities
#h: homophily factor
#delta: average density
rategen_check <- function(h, n, k, delta, N, niter, nstart) {
  classes <- rep(1:k, each = n/k)
  perm <- sample.int(n, n)
  classes <- classes[perm]
  
  W <- (1/h)*(diag(rep((h-1), k)) + outer(rep(1,k), (rep(1,k))))
  
  rb <- rbeta(n, 1, 5)
  t <- ifelse(rb < 0.05, 0.05, rb)
  
  tm <- mean(t)
  W <- (delta/(sum(W)*tm^2*(1/k^2)))*W
  
  networks <- list(0)
  for (i in 1:N) networks[[i]] <- DCBM.fast(n, k, W, t, classes)
  
  result <- ratecal(networks, h, k, classes, niter, nstart)
  
  result
}

#Calculate errors from predictive assignment using RSC for different values of p = (log m)/(log n)
ratecal <- function(networks, h, k, classes, niter, nstart) {
  N <- length(networks)
  sourceCpp("crwsample.cpp")
  
  rate1 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.80, F, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate2 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.85, F, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate3 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.90, F, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate4 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.95, F, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  
  d1 <- as.data.frame(cbind(rbind(rate1, rate2, rate3, rate4), rep("SRS", N*4), rep("RSC", N*4), rep(h, N*4), rep(c(0.80, 0.85, 0.90, 0.95), each = N)))
  names(d1) <- c("rate1", "rate2", "rate", "t1", "t2", "t3", "time", "method", "comdet", "h", "p")
  
  rate1 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.80, T, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate2 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.85, T, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate3 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.90, T, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  rate4 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- effclustering_DCBM(networks[[i]], k, 0.95, T, "np", niter, nstart))[3]
    err <- calc_error_eff(out, classes)
    c(err, out$times, unname(time))
  }))
  
  d2 <- as.data.frame(cbind(rbind(rate1, rate2, rate3, rate4), rep("RWS", N*4), rep("RSC", N*4), rep(h, N*4), rep(c(0.80, 0.85, 0.90, 0.95), each = N)))
  names(d2) <- c("rate1", "rate2", "rate", "t1", "t2", "t3", "time", "method", "comdet", "h", "p")
  
  rate1 <- t(sapply(1:N, function(i) {
    time <- system.time(out <- specclustering_DCBM(networks[[i]], k, niter, nstart))[3]
    err <- calc_error_spec(out, classes)
    c(err, 0, err, 0, sum(out$times), 0, unname(time), "full", "RSC")
  }))
  
  d3 <- as.data.frame(cbind(rbind(rate1), rep(h,N), rep(1, N)))
  names(d3) <- c("rate1", "rate2", "rate", "t1", "t2", "t3", "time", "method", "comdet", "h", "p")
  
  rbind(d1, d2, d3)
}

out <- rategen_check(h = 5, n = 100000, k = 20, delta = 0.01, N = 4, niter = 1000, nstart = 2000)
out
