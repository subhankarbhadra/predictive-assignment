function [idx,U,time_total]=consistent_wrapper2_faster(nbsAll, ClustAll,k,verbose)
T=size(nbsAll,1);
deg=sum(nbsAll');
n=size(nbsAll,2);
[~,sorted]=sort(deg,'descend');
%sorted=1:size(nbsAll,2);
iter=1;
i=1;
%Check if this has only one cluster type.
%sorted=[];
ind=1;
% for i=1:length(sorted1)
%     a=ClustAll(sorted1(i),:);
%     if(length(unique(a(a>0)))>1)
%         sorted=[sorted sorted1(i)];
%         ind=ind+1;
%     end
%     %unique(a(a>0))
%     %length(unique(a(a>0)))
%     if (ind>20)
%         break
%     end
%    
% end
%  keyboard
if nargin==3
    verbose=0;
end
adj_subsets=build_adj(nbsAll,0.9);
%[Y,N,parent,visited,score,U(iter,:)]=consistent_sequentially_patching(nbsAll,ClustAll,sorted(1),[],[],k,[],verbose);

%y=U(iter,:);

%min(max([ceil(T^(2/3)),20]),size(nbsAll,1))

%keyboard
%while iter<size(N,1)-1 && (full(sum(sum(Y(visited==1,:)>0)>0)/size(Y,2))<.9|| median(score(visited==1))<.5)
time_per_seed=0;time_align=0;
%
parfor i=1:min(max([ceil(T^(2/3)),20]),T)
    [path]=find_lenT_path(adj_subsets,sorted(i),[],[]);
    [score,all_hist,visited]=consistent_sequentially_patching_faster1(nbsAll,ClustAll,path,k,verbose);
    %[Y,N,parent,visited,score,~,all_hist]=consistent_sequentially_patching_faster(nbsAll,ClustAll,sorted(i),[],[],k,[],verbose,adj_subsets);
    y=zeros(1,n);
    a = sum(nbsAll(visited==1,:),1);
    tau=0;
    %tau=quantile(a,.01);
    %tau=max(full(mean(a)-2*std(a)),0);
    %y(a>tau)=myround(full(Y),visited,a,tau,k);
    y(a>tau)=myround_wo_hist(all_hist,a,tau,k);
    %keyboard
    scoreV(i)=median(score(visited==1));
    U(i,:)=y;
    %keyboard
end

while iter<size(U,1) %&& (sum(y>0)/length(y)<.9|| median(score(visited==1))<.7)
    %sum(visited)/numseeds
    %mean(score(visited==1))
% while iter<T
 %   i=i+1;
    iter=iter+1;
    tic;
   % [Y,N,parent,visited,score,U(iter,:)]=consistent_sequentially_patching(nbsAll,ClustAll,sorted(i),[],[],k,[],verbose);
    time_per_seed=time_per_seed+toc;
   % scoreV(iter)=median(score(visited==1));
    tic
    ap=U(iter-1,:);a=U(iter,:);
    overlap=find(ap>0&a>0);
    conf=myconfusionmat(a(overlap),ap(overlap),k,n);
    [perm]=bipartite_matching2(conf);
   % keyboard
    a(a>0)=perm(a(a>0));
    U(iter,:)=a;
    time_align=time_align+toc;
    %keyboard
   % imagesc(U)
   % pause
end
y=round(quantile(U(scoreV>=median(scoreV),:),.6));
if sum(isnan(y))>0
   y= zeros(1,length(y));
   warning('GALE Returning all zeros, since nan was detected. Check your input.');
end
%y=round(quantile(U(scoreV>quantile(scoreV,.7),:),.6));
time_total=time_align+time_per_seed/(iter-1);
%[~,index]=max(scoreV);
%idx=U(index,:);
%%% want a column vector for idx and also labels starting from 1 not from
%%% 0, hence the following two lines --- soumendu
idx = y';
if(min(idx)==0)
    idx = idx + 1;
end
%keyboard
%a=ones(1,size(A,1));
%S=find(visited==1&sum(N')>tau);
%idx=kmeans((sum(Y(visited==1,:))./sum(N(visited==1,:)))',k);