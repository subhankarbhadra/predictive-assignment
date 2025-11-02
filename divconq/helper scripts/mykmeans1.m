function [idx,C,iter] = mykmeans1(v,k,init,x)
% Uses matlab's kmeans only with several random restarts to avoid returning
% a bad solution 
% v is the vector you want to cluster based on, k is the number of clusters
% str1 and str2 are the input to the kmeans function provided by matlab

iter=1;maxiter=20; % Used for JRSS and possibly AOS (not sure about AOS)
iter=1;maxiter=20; %changing maxiter for NIPS 2015.


while (iter<=maxiter)
    if nargin==2
        %keyboard
        [idx1(iter,:),C1{iter},tmp]=kmeans(v(:,1:min(k,size(v,2))),k,'EmptyAction','drop');
    else
        [idx1(iter,:),C1{iter},tmp]=kmeans(v(:,1:min(k,size(v,2))),k,'EmptyAction','drop','Start',init');
    end
    %distVec(iter)=max(tmp);
    distVec(iter)=sum(tmp);
%     dvec(iter) = max(tmp);
%      if abs(max(tmp)-1.6648)<1e-3
%          sum(idx1(iter,:)==1)
%          keyboard
%      end
%    counts = histc(idx1(iter,:),1:k);
%     plot(idx1(iter,:),'o')
%     pause
    %dvec(iter)=min(counts);
    
    iter=iter+1;
end
[~,ind]=min(distVec);
idx = idx1(ind,:);
C=C1{ind};

%distVec
% if sum(idx==1)<.1*size(v,1) | sum(idx==2)<.1*size(v,1)
%    keyboard 
% end
idx=idx';
end


% function q=cluster_quality(A,idx)
%    v=unique(idx);
%    for i=1:length(v)
%        
%    end
% end
