function [perm]=bipartite_matching2(conf)
k=size(conf,1);
perm=zeros(k,1);
conf1=conf;

for i=1:k
    [~,ind] = max(conf1(:));
    [m,n] = ind2sub(size(conf1),ind);
    %[m n]
    perm(m)=n;
    conf1(m,:)=-1;conf1(:,n)=-1;
end