%{
Example usage: 
[A, true_id, P, Z] = cbm_parallel(8, 1000, 0.025, [0.5 0.2; 0.2 0.4], [0.7 0.3]);
%}
function [A, comm] = cbm_parallel(npar, n, rho, B, pi)

% [A, comm, P] = cbm(npar, n, rho, B, pi, true)
% [A, comm, P, Z] = cbm(npar, n, rho, B, pi, true)

if(isempty(gcp('nocreate')))
    parpool(npar);
end

k = size(B, 1);

if nargin == 4
    pi = ones(1, k) / k;
end

% if nargin < 6
%     true = 1; 
% end

% Z = sparse(n, k);
comm = randsample(k, n, true, pi);

% comm_size = floor(n * pi);
% tot = n - sum(comm_size);
% ind = randsample(k, tot);
% comm_size(ind) = comm_size(ind) + 1;
% 
% c = zeros(k + 1, 1);
% c(2:(k+1)) = cumsum(comm_size);
% 
% comm = zeros(n, 1);
% 
% for i = 1:k
%     comm((c(i) + 1):c(i+1)) = i * ones(comm_size(i));
% end

% P = rho * Z * B * Z';
% 
% if true
%     P = P - diag(diag(P));
% end

%A = sparse(n, n);
% no_edges_avg = n^2 * rho * pi * B * pi' / 2
% index = zeros(floor(1.1 * no_edges_avg), 2);
% k = 0;

% nsets = 10000;
% setsize = n / nsets;
% 
% matind = cell(nsets,1);
% 
% parfor s = 1:nsets
%     for j = ((s - 1)*setsize + 1) : (s * setsize)
%         if(j < n)
%             for i = (j + 1) : n
%                 if(rand(1) < rho * B(comm(i), comm(j)))
%                     matind{s} = [matind{s}; i, j];
%                 end
%             end
%         end
%     end
% end

matind = cell(n - 1, 1); 

parfor j = 1:n-1
    for i = (j + 1) : n
        if(rand(1) < rho * B(comm(i), comm(j)))
            matind{j} = [matind{j}; i, j];
        end
    end
end

index = vertcat(matind{:});
        
% parfor i = 2:n
%     for j = 1:(i - 1)
%         if(rand(1) < rho * B(comm(i), comm(j)))
%             A(i, j) = 1;
%         end
%     end
% end

% for i = 2:n
%     for j = 1:(i - 1)
%         if(rand(1) < rho * B(comm(i), comm(j)))
%             k = k + 1;
%             index(k, :) = [i j]; 
%         end
%     end
% end
% index = index(1:k, :);
% fprintf("Total number of edges = %d\n", k);
% 
A = sparse(index(:, 1), index(:, 2), ones(size(index, 1), 1), n, n);

A = A + A';

% toc;
end
