function [idx,U,time_total]=consistent_wrapper1(nbsAll, ClustAll,k,verbose)

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
[Y,N,parent,visited,score,U(iter,:)]=consistent_sequentially_patching(nbsAll,ClustAll,sorted(1),[],[],k,[],verbose);
scoreV(iter)=mean(score(visited==1));
y=U(iter,:);
T=size(nbsAll,1);
%min(max([ceil(T^(2/3)),20]),size(nbsAll,1))
%U(iter,:)=y;
%keyboard
%while iter<size(N,1)-1 && (full(sum(sum(Y(visited==1,:)>0)>0)/size(Y,2))<.9|| median(score(visited==1))<.5)
time_per_seed=0;time_align=0;
T1=min(max([ceil(T^(2/3)),20]),T);
parfor i=1:T1
    [Y,N,parent,visited,score,U(i,:)]=consistent_sequentially_patching(nbsAll,ClustAll,sorted(i),[],[],k,[],verbose);
    scoreV(i)=median(score(visited==1));
end

while iter<T1% && (sum(y>0)/length(y)<.9|| median(score(visited==1))<.7)
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