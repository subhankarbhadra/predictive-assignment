function B = rproj(A, d)
    n = size(A, 2);
    B = A*randn(n,d)/sqrt(d);
end
