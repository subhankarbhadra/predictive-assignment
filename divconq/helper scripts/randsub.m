function [N, ind] = randsub(A, m)
n = size(A, 1);

ind = randsample(n, m);
N = A(ind, ind);
end
