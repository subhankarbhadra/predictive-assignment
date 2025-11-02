function [path,visited]=find_lenT_path(adj_graph,parent,path,visited)
T=size(adj_graph,1);


thresh=.55;

%first=N{1};t=sparse(n,1);t(first)=1;
if isempty(path)
    visited=zeros(1,T);
    path=[];
end



%for i=1:num_seeds
%visited(parent)
if visited(parent)>0
    %disp('blah')
    return;
end

visited(parent)=1;
path=[path parent];
%[~,inds]=sort(N*N(parent,:)','descend');
%nbs1=inds;
%ov=N*N(parent,:)';
%median(ov)
%nbs1=find(squeeze(ov>quantile(ov,.9)));
nbs1=find(adj_graph(parent,:));
%quantile(ov,.99)
%sum(visited(nbs)==0)
%sum(visited(nbs1)>0)/length(nbs1)
%parent
%y=zeros(1,n);


    
%tau=quantile(visited_stat,.01);

for j=nbs1
    
    if visited(j)>0||visited(j)<-5
        continue
    end
    
    [path,visited]=find_lenT_path(adj_graph,j,path,visited);

    
end

%end
%plot(sum(Y(visited==1,:))./sum(N(visited==1,:)));
%keyboard
end