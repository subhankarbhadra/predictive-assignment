function adj_subsets=build_adj(nbsAll,p)
T=size(nbsAll,1);
NN=nbsAll*nbsAll';
adj_subsets=zeros(T,T);
parfor i=1:T
    adj_subsets(i,:)=NN(i,:)>quantile(NN(i,:),p,2);
end
%nbs1=find(squeeze(ov>quantile(ov,.9)));