%%%%%%%%%%
% Authors : Soumendu Sundar Mukherjee, Purnamrita Sarkar
% License: LGPLv3
% Description: Wrapper script for PACE and GALE
%%%%%%%%%% 
%
% A = input adjacency matrix
%
% COV = node covariates (optional), for not specifying use [].
%
% T = number of subgraphs used
%
% method = 'PACE' or 'GALE', default is PACE
%
% recovery_method = 'naive' or 'kmeans' or 'spectral' or 'na', used only with PACE, for
% GALE use 'na', default is 'naive'
%
% subgraph_type = 'random' or 'nbhood' or 'onion', default is 'nbhood'
%
% size = size of subgraph, if method == 'random', else size = # of hops to
% be taken for nbhood
% 
% no_of_blocks = self explanatory, defaults to 2 if not specified
%
% algorithm = 'sp' for spectral clustering, 'rsp' for regularised laplacian spectral clustering, 'asp'
% for adjacency spectral clustering, 'sdp' for semi definite programming,
% 'pl' for profile likelihood, 'psl' for pseudo likelihood, 'mfl' for
% mean-field likelihood, 'custom' for user supplied method (pass a function handle in this case)
%
% optargs = specify other optional arguments as a struct object: optargs.tau_method = 'absolute' or 'dynamic', optargs.tau_param, optargs.user_init = true of false, optargs.start, optargs.row_normalization, optargs.true_id ...
%
%
% Example 1:
% [C, Z, sigma] = divconq_parallel(A, [], 100, 'onion', 3, 'PACE', 'kmeans', 2, 'asp',optargs)
%
%
%

function [comm, Z, C, avg_subgraph_size, timecons,toStore] = divconq_parallel(A, COV, T, subgraph_type, subgraph_param, method, recovery_method, no_of_blocks, algorithm, optargs)
n = size(A, 1);     % size of the network
%n = length(A);
K = no_of_blocks;
name=sprintf('n_%d_rho_%.3f_bal_%d_T%d_%s_%d_%s_%s_k%d_alg_%s.mat',n,optargs.rho,optargs.bal,T,subgraph_type, subgraph_param,method, recovery_method,no_of_blocks,algorithm);
%C = zeros(n, n);
Z_both = struct;
Z_both.PACE = zeros(n, K);
Z_both.GALE = zeros(n, K);
comm_both = struct;
comm_both.PACE = zeros(n, 1);
comm_both.GALE = zeros(n, 1);
%sigma = zeros(n,1);

mstar = 5;

nbsAll = zeros(T, n);
ClustAll = zeros(T, n);
subgraph_size = zeros(T,1);
roots = randsample(n, T);

switch(subgraph_type)
    case('random')
        sampler = @(x,y,z)randsub(x,z);
    case('nbhood')
        sampler = @(x,y,z)hhop(x,y,z);
    case('onion')
        sampler = @(x,y,z)onion(x,y,z);
end

c = 0;

switch(algorithm)
    case('sp')
        row_normalization = optargs.row_normalization;
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [~, comm_sub] = spectral(As, K, 'unregLaplacian', row_normalization);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('rsp')
        row_normalization = optargs.row_normalization;
        reg_type = optargs.reg_type1;
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [~, comm_sub] = spectral(As, K, 'regLaplacian', row_normalization, reg_type);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('asp')
        row_normalization = optargs.row_normalization;
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [~, comm_sub] = spectral(As, K, 'adjacency', row_normalization);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('sdp')

        opts.T=300;opts.tol=1e-3;opts.rho=1;
        opts.report_interval=1;
        opts.quiet=1;
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            %length(verts)
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [X,~] = admm_imb(As,K,1/(size(As,1)/2),opts);
                [~, comm_sub] = spectral(X, K, 'regLaplacian', true);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('plh')          % Profile likelihood (heuristic)
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [~, comm_sub, ~] = heuristic_max_profile_likelihood(As,K);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('pl')          % Profile likelihood (tabu search)
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [~, comm_sub, ~] = tabuDiffStart(As,K);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('psl')         % Pseudo-likelihood
        true_id = optargs.true_id;
        % options for the init method and cpl/upl
        init_opts = struct('verbose',false);
        Tpsl = 20;
        cpl_opts = struct('verbose',false,'delta_max',0.1,'itr_num',Tpsl,'em_max',80,'track_err',false);
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                [e, ~] = initLabel5b(As, K, 'scp', init_opts);
                [comm_sub, ~, ~, ~] = cpl4c(As, K, e, true_id, 'cpl', cpl_opts);
                if(length(verts)~=length(comm_sub))
                    fprintf('Zero degree nodes detected. Remove these before using divconq.');
                    keyboard
                end
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;
    case('mfl')
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                if optargs.user_init
                    start = optargs.start;
                    start1 = struct;
                    start1.Z = start.Z(verts,:);
                    start1.B = start.B;
                    [Zs, ~] = meanfield_generic(As, K, 1, start1);
                else
                    [Zs, ~] = meanfield_generic(As, K, 1);
                end
                [~, comm_sub] = max(Zs, [], 2);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('mflcov')         % Mean-field likelihood with covariates
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                if optargs.user_init
                    start = optargs.start;
                    maxiter1 = start.maxiter1;
                    maxiter2 = start.maxiter2;
                    stepsize = start.stepsize;
                    start2 = struct;
                    start2.psi = start.psi(verts,:);
                    start2.theta = start.theta;
                    start2.beta = start.beta;
                    start2.pi = start.pi;
                    [Zs, ~] = meanfield_cov(As,COV(verts,:),K,start2,maxiter1,maxiter2,stepsize);
                else
                    [Zs, ~] = meanfield_generic(As, K, 1);
                end
                [~, comm_sub] = max(Zs, [], 2);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;

    case('custom') 
        tstart = tic;
        parfor t = 1:T
            [As, verts] = sampler(A, roots(t), subgraph_param);
            % remove zero degree nodes
            zdNodes = sum(As,2) == 0;
            As = As(~zdNodes, ~zdNodes);
            verts = verts(~zdNodes);
            if(length(verts) >= mstar)
                c = c+1;
                subgraph_size(t) = length(verts);
                X = sdpcvx(As, K, 'SDP-3');
                [~, comm_sub] = spectral(X, K, 'regLaplacian', true);
                v1 = zeros(1,n);
                v1(verts) = 1;
                nbsAll(t, :) = v1;

                v2 = zeros(1, n);
                v2(verts) = comm_sub';
                ClustAll(t, :) = v2;
            end 
        end
        tot_time = toc(tstart);
        timecons.GALE_subgraphs_total = tot_time;
        timecons.GALE_subgraphs_avg = tot_time/T;
    otherwise
        fprintf('Error: algorithm not supported!\n');
end

avg_subgraph_size = sum(subgraph_size)/c;
fprintf('\n Average subgraph size = %3.2f\n', avg_subgraph_size);
N = zeros(n, n);
Ctemp = N;
C = N;

percflag = true;
if strcmp(optargs.tau_method, 'absolute')
    tau = optargs.tau_param;
    percflag = false;
elseif strcmp(optargs.tau_method, 'dynamic')
    perc = optargs.tau_param;
else
    perc = 0.4;
end

for t = 1:T
    verts = find(nbsAll(t,:) == 1);
    comm = ClustAll(t, verts);
    Zs = zeros(length(verts),K); 
    for i = 1:length(verts)
        Zs(i,comm(i)) = 1;
    end
    Zext = zeros(n, K);
    Zext(verts, :) = Zs;
    N(verts, verts) = N(verts, verts) + 1;
    %W(verts, verts) = W(verts, verts) + length(verts);
    Ctemp = Ctemp + Zext * Zext';
    %CWtemp = CWtemp + length(verts) * Zext * Zext';
end


tstart1 = tic;
if percflag
    tau = quantile(N(N > 0), perc);
end
%         C = N; % hack
%         C(C >= tau) = Ctemp(N >= tau) ./ N(N >= tau);
for i = 1:n
    for j = 1:n
        if(N(i,j) >= tau)
            C(i,j) = Ctemp(i,j)/N(i,j);
        end
%         if(W(i,j) >= Wtau)
%             CW(i,j) = CWtemp(i,j)/W(i,j);
%         end
    end
end
timecons.PACE_thresholding_C = toc(tstart1);
%keyboard
% Recoveing Z from C
tstart2 = tic;
if(strcmp(recovery_method, 'naive'))
    C_proj = rproj(C, 50*floor(log(n))); % random projection
    comm_both.PACE = naivecluster(C, C_proj, K); % naive clustering
elseif(strcmp(recovery_method, 'kmeans'))
    C_proj = rproj(C, 50*floor(log(n))); % random projection
    comm_both.PACE = mykmeans1(C_proj, K);   % kmeans clustering
    %keyboard
else
    [~, comm_both.PACE] = spectral(C, K, 'unregLaplacian', 'true');
    end

for i = 1:n
    Z_both.PACE(i, comm_both.PACE(i)) = 1;
end
timecons.PACE_recovering_Z = toc(tstart2);


tstart3 = tic;
[comm_both.GALE, ~] = consistent_wrapper1(nbsAll, ClustAll, K);
timecons.GALE_patching_step = toc(tstart3);
for i = 1:n
    Z_both.GALE(i, comm_both.GALE(i)) = 1;
end


switch(method)
    case('PACE')
        comm = comm_both.PACE;
        Z = Z_both.PACE;        
    case('GALE')
        comm = comm_both.GALE;
        Z = Z_both.GALE;                
        C = [];
    otherwise
        comm = comm_both;
        Z = Z_both;
end
%keyboard

toStore.nbsAll=nbsAll;
toStore.ClustAll=ClustAll;
toStore.C=C;
toStore.N=N;
toStore.comm_both=comm_both;
toStore.time=timecons;
toStore.truth=optargs.true_id;
timecons;
if ~exist(name, 'file')==1
    disp('saving files')
    save(name,'-struct','toStore')
end
%save name.mat toStore
end