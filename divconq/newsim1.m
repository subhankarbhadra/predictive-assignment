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

%%
rng(12345);
nworkers = 18;
N = 50;
rates = zeros(N, 8);
times = zeros(N, 8);

for i = 1:N
    [A, comm] = cbm_parallel(nworkers, n, rho, B, pi);
    n = size(A, 1); K = max(comm);

    time = tic;
    [~, comm_est] = spectral(A, K, 'unregLaplacian', 'true');
    t1 = toc(time);
    err1 = cluster_acc(comm, comm_est);

    time = tic;
    [~, comm_est] = spectral(A, K, 'regLaplacian', 'true', 'Amini');
    t2 = toc(time);
    err2 = cluster_acc(comm, comm_est);

    [pred,~,~,~,tcon] = divconq_parallel(A, [], 100, 'random', 1500, 'both', 'spectral', K, 'sp', ...
        struct('rho', 0.002, 'bal', true, 'row_normalization', 'true', 'tau_method', 'custom', 'true_id', comm));

    err3_p = cluster_acc(comm, pred.PACE);
    err3_g = cluster_acc(comm, pred.GALE);

    t3_p = tcon.GALE_subgraphs_total + tcon.PACE_thresholding_C + tcon.PACE_recovering_Z;
    t3_g = tcon.GALE_subgraphs_total + tcon.GALE_patching_step;

    [pred,~,~,~,tcon] = divconq_parallel(A, [], 100, 'random', 1500, 'both', 'spectral', K, 'rsp', ...
        struct('rho', 0.002, 'bal', true, 'row_normalization', 'true', 'tau_method', 'custom', 'true_id', comm, ...
        'reg_type1', 'Amini'));

    err4_p = cluster_acc(comm, pred.PACE);
    err4_g = cluster_acc(comm, pred.GALE);

    t4_p = tcon.GALE_subgraphs_total + tcon.PACE_thresholding_C + tcon.PACE_recovering_Z;
    t4_g = tcon.GALE_subgraphs_total + tcon.GALE_patching_step;

    time = tic;[pred2, ~] = effclustering(A, K, 0.5, 'sp', 'cc');t5 = toc(time);
    err5 = cluster_acc(comm, pred2);

    time = tic;[pred2, ~] = effclustering(A, K, 0.5, 'rsp', 'cc', 'Amini');t6 = toc(time);
    err6 = cluster_acc(comm, pred2);

    rates(i, :) = [err1 err2 err3_p err3_g err4_p err4_g err5 err6];
    times(i, :) = [t1 t2 t3_p t3_p t4_g t4_g t5 t6];
end

%%
% Average community detection error and runtime for SC, RSC-A, SC + PACE,
% SC + GALE, RSC-A + PACE, RSC-A + GALE, SC + Predictive Assignment, and 
% RSC-A + Predictive Assignment 
mean(rates, 1)
std(rates, 1)
mean(times, 1)
