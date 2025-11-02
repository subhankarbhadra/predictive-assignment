%{
Arguments:
    npar: number of parallel cores to use
    B: block probability matrix
    comm: community membership vector
%}
function [A] = cbm_parallel_v2(npar, B, comm)

if(isempty(gcp('nocreate')))
    parpool(npar);
end

n = length(comm);

matind = cell(n - 1, 1); 

parfor j = 1:(n - 1)
    for i = (j + 1) : n
        if(rand(1) < B(comm(i), comm(j)))
            matind{j} = [matind{j}; i, j];
        end
    end
end

index = vertcat(matind{:});
        
A = sparse(index(:, 1), index(:, 2), ones(size(index, 1), 1), n, n);

A = A + A';

end
