%%%%%%%%%%
% Authors : Soumendu Sundar Mukherjee, Purnamrita Sarkar
% License: LGPLv3 (inclusive of all the scripts in the 'helper scripts' folder)
% Description: Large scale simulations
%%%%%%%%%% 

% Note: Keep 'patching_large_sim.m' and 'helper scripts' in the same directory

% Load helper scripts
addpath ./'helper scripts'/

% Set seed for reproducibility
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

% Set up parallel pool
npar = 100;

pool = gcp('nocreate');

if isempty(pool)
	parpool(npar);
elseif pool.NumWorkers ~= npar
	delete(pool);
	parpool(npar);
end

% Simulation parameters
n = 10^7;

m = floor(n^(2/3));

T = 10 * floor(n / m);

K = 20;

comm = ones(n / K, 1);
for k = 2:K
    comm = [comm; k * ones(n / K, 1)];
end

a = 500;
b = 10;

p = a / m;
q = b / m;

B = (p - q) * eye(K) + q * ones(K, K);

deg_sub = (a + b * (K - 1)) / K;
deg = (n / m) * deg_sub;

fprintf('\nAverage degree per subgraph = %.2f', deg_sub);
fprintf('\nAverage degree = %.2f\n', deg);

% Sample subgraphs, while taking care of overlaps
subgraphs = cell(T, 1);
vertices = cell(T, 1);

fprintf('Generating subgraph %d...\n', 1);
tic;
vertices{1} = randsample(n, m);
subgraphs{1} = cbm_parallel_v2(npar, B, comm(vertices{1}));
toc;

for t = 2:T
    fprintf('\nGenerating subgraph %d...\n', t);
    tic;
    vertices{t} = randsample(n, m);
    subgraphs{t} = cbm_parallel_v2(npar, B, comm(vertices{t})); % wasteful
    toc;
    remaining = vertices{t};
    fprintf('Taking care of overlaps...\n');
    tic;
    for s = 1:(t - 1)
        overlap = intersect(vertices{s}, remaining);
        ind1 = find(ismember(vertices{s}, overlap) == 1);
        ind2 = find(ismember(vertices{t}, overlap) == 1);
        subgraphs{t}(ind2, ind2) = subgraphs{s}(ind1, ind1);
        remaining = setdiff(remaining, overlap);
    end
    toc;
end

covered = unique( vertcat(vertices{:}) );
fprintf('\nFraction of sampled vertices = %.2f\n', length(covered)/n);

nbsAll = sparse(T, n);
ClustAll = sparse(T, n);

mstar = 10;
thresh = 1;

c = 0; % accepted subgraphs out of T
avg_subgraph_size = 0;


err_sub = zeros(T, 1);

tstart = tic;

parfor t = 1:T
    As = subgraphs{t};
    verts = vertices{t}
    [As, ~, ~, I] = process_real_graph(As, As, As, thresh); % largest connected component
	verts = verts(I);
    if(length(I) >= mstar)
        c = c+1;
        avg_subgraph_size = avg_subgraph_size + length(verts);
        [~, comm_sub] = spectral(As, K, 'regLaplacian', 1, 'Rohe');
        error_sub = cluster_acc(comm_sub, comm(verts));
        err_sub(t) = error_sub;
        v1 = zeros(1, n);
        v1(verts) = 1;
        nbsAll(t, :) = v1;
        v2 = zeros(1, n);
        v2(verts) = comm_sub;
        ClustAll(t, :) = v2;
    end
end

tot_time = toc(tstart);

timetaken.subgraphs_total = tot_time;
timetaken.subgraphs_avg = tot_time/T;

avg_subgraph_size = avg_subgraph_size/c;
fprintf('Average subgraph size = %.f\n', avg_subgraph_size);

%%%% Save workspace %%%%
filename_1 = strcat('m', num2str(m), '-K', num2str(K), '-T', num2str(T), '.mat');
save(filename_1, 'npar', 'm', 'T', 'n', 'K', 'B', 'ClustAll', 'nbsAll', 'avg_subgraph_size', 'timetaken', 'err_sub', 'comm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(gcp('nocreate'))
	parpool(npar);
end

%%%% GALE %%%%
fprintf('Performing GALE specific computations...\n');
tstart1 = tic;
[comm_both.GALE, ~] = consistent_wrapper2_faster(nbsAll, ClustAll, K, 0);
timetaken.GALE_patching_step = toc(tstart1);

err_GALE = cluster_acc(comm_both.GALE, comm);

%%%% Save workspace %%%%
filename_2 = strcat('m', num2str(m), '-K', num2str(K), '-T', num2str(T), '-complete-', num2str(npar), 'cores.mat');
save(filename_2, 'npar', 'm', 'T', 'n', 'K', 'B', 'ClustAll', 'nbsAll', 'timetaken', 'err_sub', 'err_GALE', 'comm', 'comm_both');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stdout
fprintf('\n-------------------------------------------------------------\n');
fprintf('   Algorithm  |  Misclustering error  |  Time taken (seconds)\n');
fprintf('-------------------------------------------------------------\n');
fprintf('        GALE  |         %.4f        |        %.2f\n', err_GALE, (timetaken.subgraphs_total + timetaken.GALE_patching_step));
fprintf('-------------------------------------------------------------\n\n');
