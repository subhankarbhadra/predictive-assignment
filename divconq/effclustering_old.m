function [comm_pred, time] = effclustering_old(A, K, p, sampling, comdet, row_normalization, method, reg_type)
n = size(A, 1); % size of the network
m = floor(n*p);

switch(sampling)
    case('SRS')
        sampler = @(x,y,z)randsample(y, z);
    case('RWS')
        fprintf('Error: algorithm not supported!\n');
end
subnet = sampler(A, n, m);
As = A(subnet, subnet);
% remove zero degree nodes
zdNodes = sum(As,2) == 0;
As = As(~zdNodes, ~zdNodes);
subnet = subnet(~zdNodes);
rest = setdiff(1:n, subnet);
comm_pred = zeros(n, 1);

switch(comdet)
    case('sp')
        tstart = tic;
        [~, comm_sub] = spectral(As, K, 'unregLaplacian', row_normalization);
        tot_time = toc(tstart);
        time.comdet = tot_time;

    case('rsp')
        tstart = tic;
        [~, comm_sub] = spectral(As, K, 'regLaplacian', row_normalization, reg_type);
        tot_time = toc(tstart);
        time.comdet = tot_time;

    case('basp')
        tstart = tic;
        D = diag(sum(A, 1));
        Bs = A(subnet, :)*A(:, subnet) - D(subnet, subnet);
        [~, comm_sub] = spectral(Bs, K, 'unregLaplacian', row_normalization);
        tot_time = toc(tstart);
        time.comdet = tot_time;
        fprintf('Error: algorithm not supported!\n');
end

comm_pred(subnet) = comm_sub;
mr = histc(comm_sub(:),unique(comm_sub));
switch(method)
    case('cc')
        tstart2 = tic;
        theta = zeros(n, K);
        for k = 1:K
            theta(:, k) = sum(A(:, subnet(comm_sub == k)), 2)/mr(k); 
        end
        a = sum(theta.^2, 1);
        [i,j] = find(A);
        func = @(x) {x};
        nbr = splitapply(func, vertcat(i, [1:n]'), vertcat(j, [1:n]'));
        %for i = rest
        %    [~,comm_pred(i)] = min(-2*sum(A(:, i).*theta, 1) + a);
        %end
        for i = rest
            [~,comm_pred(i)] = min(-2*sum(theta(nbr{i},:), 1) + a);
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