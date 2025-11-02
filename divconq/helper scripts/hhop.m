% returns the h-hop neighbourhood of vertex root
function [N, nbs] = hhop(A, root, h)
n = size(A, 1);

Ah = speye(n,n);
Bh = Ah;
for j = 1:h
    Ah = Ah*A;
    Bh = Bh + Ah;
end

nbs = find(Bh(root,:));
N = A(nbs, nbs);
end