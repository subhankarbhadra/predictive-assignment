source("predinf_sbm.R")

h = 3; n = 50000; k = 15; delta = 0.01
# h = 3; n = 100000; k = 20; delta = 0.01
# h = 3; n = 150000; k = 20; delta = 0.01
# h = 3; n = 200000; k = 20; delta = 0.01

classes <- sort(rep(1:k, ceiling(n/k))[1:n])
perm <- sample.int(n, n)
classes <- classes[perm]

W <- (delta*k/(h + k - 1))*(diag(rep((h-1), k)) + outer(rep(1,k), (rep(1,k))))

A <- as_adjacency_matrix(sample_sbm(n, W, as.numeric(table(classes))))[perm, perm]
saveRDS(A, "sbmnet_n50k_k15_h3.rds")

# Open a clean R session and run the following code
#
source("predinf_sbm.R")

A <- readRDS("sbmnet_n50k_k15_h3.rds")
k <- 15; niter <- 1000; nstart <- 2000

peakRAM(effclustering(A, k, 0.70, "cc", T, niter, nstart))
peakRAM(effclustering(A, k, 0.75, "cc", T, niter, nstart))
peakRAM(effclustering(A, k, 0.80, "cc", T, niter, nstart))
peakRAM(effclustering(A, k, 0.85, "cc", T, niter, nstart))
peakRAM(specclustering(A, k, T, niter, nstart))
