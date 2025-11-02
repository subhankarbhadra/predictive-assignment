% Clustering based on predictive assignment (closest community approach)
%
% Arguments:
%  A: n x n symmetric binary adjacency matrix
%  K: Integer, number of communities
%  p: Numeric, (m/n) where m is the size of the subgraph
%  comdet: Community detection method on subgraph 
%          ("sc" for spectral clustering, 
%           "rsc" for regularized spectral clustering, 
%           "basc" for bias-adjusted spectral clustering)
%  method: Predictive assignment approach 
%          ("cc" for closest community approach)
%  reg_type: Type of regularization 
%            ("rohe" for Rohe-type regularization,
%             "amini" for Amini-type regularization)
%
% Returns: A list
%  comm_pred: Integer vector, estimated community assignments
%  time: Numeric, total runtime
%
function [comm_pred, time] = effclustering(A, K, p, comdet, method, reg_type)
n = size(A, 1); % size of the network
m = floor(n*p);

sampler = @(x,y,z)randsample(y, z);
subnet = sampler(A, n, m);
As = A(subnet, subnet);
% remove zero degree nodes
zdNodes = sum(As,2) == 0;
As = As(~zdNodes, ~zdNodes);
subnet = subnet(~zdNodes);
rest = setdiff(1:n, subnet);
m = length(subnet);
comm_pred = zeros(n, 1);

switch(comdet)
    case('sp')
        tstart = tic;
        [~, comm_sub] = spectral(As, K, 'unregLaplacian', 'true');
        tot_time = toc(tstart);
        time.comdet = tot_time;

    case('rsp')
        tstart = tic;
        [~, comm_sub] = spectral(As, K, 'regLaplacian', 'true', reg_type);
        tot_time = toc(tstart);
        time.comdet = tot_time;

    case('basp')
        tstart = tic;
        D = diag(sum(A, 1));
        Bs = A(subnet, :)*A(:, subnet) - D(subnet, subnet);
        [~, comm_sub] = spectral(Bs, K, 'unregLaplacian', 'true');
        tot_time = toc(tstart);
        time.comdet = tot_time;
        fprintf('Error: algorithm not supported!\n');
end

comm_pred(subnet) = comm_sub;
mr = histc(comm_sub(:),unique(comm_sub));
switch(method)
    case('cc')
        tstart2 = tic;
        theta = zeros(n-m, K);
        for k = 1:K
            theta(:, k) = sum(A(rest, subnet(comm_sub == k)), 2)/mr(k); 
        end
        a = sum(theta.^2, 1);
        [i,j] = find(A(rest, rest));
        func = @(x) {x};
        nbr = splitapply(func, vertcat(i, [1:(n-m)]'), vertcat(j, [1:(n-m)]'));
        %for i = rest
        %    [~,comm_pred(i)] = min(-2*sum(A(:, i).*theta, 1) + a);
        %end
        for i = 1:(n-m)
            [~,comm_pred(rest(i))] = min(-2*sum(theta(nbr{i},:), 1) + a);
        end
        tot_time = toc(tstart2);
        time.classification = tot_time;
    case('np')
        fprintf('Error: algorithm not supported!\n');
end
%e = cluster_acc(comm_pred, true_comm);
%e1 = cluster_acc(comm_sub, true_comm(subnet));
%e2 = cluster_acc(comm_pred(rest), true_comm(rest));
time
end