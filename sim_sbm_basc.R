source("predinf_sbm.R")

# Simulate networks from SBM and apply predictive assignment
# Arguments:
#  h: Numeric, homophily factor
#  n: Integer, number of nodes
#  k: Integer, number of communities
#  delta: Numeric, average density of networks
#  N: Integer, number of replications
#  niter: Integer, maximum number of iterations for k-means
#  nstart: Integer, number of initializations for k-means
#
# Returns: 
#  Output of ratecal() for the simulated networks
#
rategen_check <- function(h, n, k, delta, N, niter, nstart) {
  classes <- sort(rep(1:k, ceiling(n/k))[1:n])
  perm <- sample.int(n, n)
  classes <- classes[perm]
  
  W <- (delta*k/(h + k - 1))*(diag(rep((h-1), k)) + outer(rep(1,k), (rep(1,k))))
  
  networks <- list(0)
  for (i in 1:N) networks[[i]] <- as_adjacency_matrix(sample_sbm(n, W, as.numeric(table(classes))))[perm, perm]
  
  result <- ratecal(networks, h, k, classes, niter, nstart)
  result
}

# Calculate errors from predictive assignment using BASC for different values of p = (log m)/(log n)
#
# Arguments:
#  networks: List of adjancency matrices
#  h: Numeric, homophily factor
#  k: Integer, number of communities
#  classes: Integer vector, true community labels
#  niter: Integer, maximum number of iterations for k-means
#  nstart: Integer, number of initializations for k-means
#
# Returns: A matrix consisting of the community detection errors and runtimes for predictive assignment using BASC for different values of p = (log m)/(log n)
#  Columns 1-3: Community detection error for subgraph nodes, remaining nodes, and all nodes 
#  Columns 4-7: Runtimes for steps 1-3 of predictive assignment and total runtime
#  Columns 8-11: Attributes used for summarizing results
#
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

# To reproduce the bias-adjusted spectral clustering results in Tables 2 and 3,
# run rategen_check() for the different choices of (h, n, k), using N=30 replications.

# For quicker demonstration, use a smaller number of replications (e.g., N=4).
# Since N=30 replications is computationally intensive, 
# we ran this script 15 times in parallel (as 15 separate jobs), each with N=2 replications, and combined results.

N = 4
h = 3; n = 50000; k = 15
# h = 3; n = 100000; k = 20
# h = 3; n = 150000; k = 20
# h = 3; n = 200000; k = 20
out <- rategen_check(h, n, k, delta = 0.01, N, niter = 1000, nstart = 2000)
for (i in c(1:7, 10, 11)) out[[i]] <- as.numeric(out[[i]])
out

# Summarizing results
# Average community detection errors for subgraph nodes, remaining nodes, and all nodes, and average runtime
out %>% group_by(comdet, h, p) %>% 
  summarise(avg.rate1 = mean(rate1)*100,
            sd.rate1 = sd(rate1)*100, avg.rate2 = mean(rate2)*100,
            sd.rate2 = sd(rate2)*100, avg.rate = mean(rate)*100,
            sd.rate = sd(rate)*100, avg.time = mean(time)) %>%
  as_tibble()

# Note: memory consumption is not evaluated here.
# Running additional processes, such as simulating networks, can affect memory computations, 
# making it difficult to compute the memory usage solely due to community detection.
# Instead, simulated networks are saved, and then loaded in a clean R session for more reliable memory computations.
# See peakram_basc_SBM.R for demonstration.