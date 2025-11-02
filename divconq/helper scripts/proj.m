function Y = proj(X, s, method)

switch(method)
    case('random')
        n = size(X, 2);
        Y = X * randn(n, s)/sqrt(s);
    case('pca')
        [V, ~] = eigs(X, s, 'largestreal');
        Y = X * V;
end
    