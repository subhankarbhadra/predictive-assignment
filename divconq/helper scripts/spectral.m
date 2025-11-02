function [Z, comm, iter] = spectral(A, K, method, row_normalized, reg_type, reg)

n = size(A,1);

if nargin == 2
    method = 'unregLaplacian';
    row_normalized=0;
end

switch(method)
    case('unregLaplacian')
        deg = A*ones(n,1);
        deginvsqrt = sqrt(deg);
        ind = (deginvsqrt > 0);
        deginvsqrt(ind) = deginvsqrt(ind).^(-1);
        Y = zeros(n, 1);
        Dinvsqrt = spdiags([Y deginvsqrt Y], -1:1, n, n);
        L = Dinvsqrt*A*Dinvsqrt;
        %L=.5*(L+L');
        %size(L)
        %keyboard
        %disp('starting svd')
        %if (size(L,1)<10)
        %    [V,E] = svd(L); % Top K eigenvalues and eigenvectors of L
        %    [~,inds]=sort(abs(diag(E)),'descend');
        %    %diag(E)
        %    V = V(:, inds(1:K));
        %    %E = E(inds(1:K), inds(1:K));%diag(E)
        %else
        %    [V, ~] = svds(L,K);
        %end
        %[V, ~] = svds(L, K, 'L');
        [V, ~] = svds(L, K, 'largest','SubspaceDimension', 100);
        %disp('finished svd')
    case('regLaplacian')
        if(nargin == 5)
            reg = mean(sum(A)); % Regularization parameter set to average degree
        end
        switch(reg_type)
            case('Amini')
                Areg = A + 0.25*reg/n; % Regularization of Amini et al. (warning: this produces a dense matrix).
                deg = Areg*ones(n, 1);
                deginvsqrt = sqrt(deg);
                ind = (deginvsqrt > 0);
                deginvsqrt(ind) = deginvsqrt(ind).^(-1);
                Y = zeros(n, 1);
                Dinvsqrt = spdiags([Y deginvsqrt Y], -1:1, n, n);
                %L = Dinvsqrt*Areg*Dinvsqrt;
                L = Dinvsqrt*A*Dinvsqrt + reg*sum(deginvsqrt)^2/n;
            case('Rohe')
                %tstart = tic;
                deg = A*ones(n, 1) + reg; % Regularization of Qin and Rohe
                deginvsqrt = deg.^(-1/2);
                Y = zeros(n, 1);
                Dinvsqrt = spdiags([Y deginvsqrt Y], -1:1, n, n);
                L = Dinvsqrt*A*Dinvsqrt;
                %tot = toc(tstart);
                %time.laplace_construct = tot
        end
                
        %L=.5*(L+L');
        %size(L)
        %keyboard
        %tstart2 = tic;
        if(size(L, 1) < 10)
            [V, E] = svds(L); % Top K eigenvalues and eigenvectors of L
            [~, inds] = sort(abs(diag(E)), 'descend');
            %diag(E)
            V = V(:, inds(1:K));
            %E = E(inds(1:K), inds(1:K));%diag(E)
        else
            %[V, ~] = svds(L, K, 'L');
            [V, ~] = svds(L, K, 'largest','SubspaceDimension', 100);
        end
        %tot2 = toc(tstart2);
        %time.SVD = tot2
    case('adjacency')
        %[V, ~] = svds(A, K, 'L');
        [V, ~] = svds(A, K, 'largest','SubspaceDimension', 100);
end
%row_normalized
%tstart3 = tic;
if(row_normalized)
    for i = 1:n
        if(norm(V(i, :)) > 0)
            V(i, :) = V(i, :)/norm(V(i, :)); % row normalization useful for degree corrected models.
        end
    end
end
%tot3 = toc(tstart3);
%time.rownorm = tot3

%tstart4 = tic;
[comm,~,iter] = mykmeans1(V, K);
%tot4 = toc(tstart4);
%time.Kmeans = tot4

Z = sparse(1:n, comm, ones(n, 1), n, K);

% Z = [];
% Z = sparse(n, K);
% for i = 1:n
%    Z(i, comm(i)) = 1;
% end
%keyboard
end
