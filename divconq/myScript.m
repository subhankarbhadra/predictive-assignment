addpath ./'helper scripts'/
%%
n = 5000;
rho = 0.002;
p = 1;
r = 0.01;
q = r*p;
K = 2;
B = (p - q) * eye(K) + q * ones(K);
pi = [.2 .8];
avg_deg = (n - 1) * rho * pi * B * pi';

nworkers = 18;

[A, comm] = cbm_parallel(nworkers, n, rho, B, pi);
n = size(A, 1); K = max(comm);

time = tic;
[~, comm_est] = spectral(A, K, 'unregLaplacian', 'true');
t1 = toc(time);
err1 = cluster_acc(comm, comm_est);
err1
