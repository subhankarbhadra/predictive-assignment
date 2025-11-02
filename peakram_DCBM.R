source("predinf_dcbm.R")

h = 3; n = 100000; k = 20; delta = 0.01
# h = 5; n = 100000; k = 20; delta = 0.01

classes <- rep(1:k, each = n/k)
perm <- sample.int(n, n)
classes <- classes[perm]

W <- (1/h)*(diag(rep((h-1), k)) + outer(rep(1,k), (rep(1,k))))

rb <- rbeta(n, 1, 5)
t <- ifelse(rb < 0.05, 0.05, rb)

tm <- mean(t)
W <- (delta/(sum(W)*tm^2*(1/k^2)))*W

A <- DCBM.fast(n, k, W, t, classes)
saveRDS(A, "dcbmnet_n100k_k20_h3.rds")

# Open a clean R session and run the following code
#
source("predinf_dcbm.R")

A <- readRDS("dcbmnet_n100k_k20_h3.rds")
k <- 20; niter <- 1000; nstart <- 2000

sourceCpp("crwsample.cpp")

peakRAM(effclustering_DCBM(A, k, 0.80, T, "np", niter, nstart))
peakRAM(effclustering_DCBM(A, k, 0.85, T, "np", niter, nstart))
peakRAM(effclustering_DCBM(A, k, 0.90, T, "np", niter, nstart))
peakRAM(effclustering_DCBM(A, k, 0.95, T, "np", niter, nstart))
peakRAM(specclustering_DCBM(A, k, niter, nstart))
