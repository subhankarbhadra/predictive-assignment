function [Y,N,parent,visited,score,y]=consistent_sequentially_patching(N,Y,parent,visited,score,k,y,verbose)
thresh=.55;
num_seeds=size(N,1);

n=size(Y,2);

if nargin==7
    verbose=0;
end
%first=N{1};t=sparse(n,1);t(first)=1;
if length(visited)==0
    visited=zeros(1,num_seeds);
    score=visited;
end



%for i=1:num_seeds
%visited(parent)
if visited(parent)>0
    %disp('blah')
    return;
end

visited(parent)=1;
%[~,inds]=sort(N*N(parent,:)','descend');
%nbs1=inds;
ov=N*N(parent,:)';
%median(ov)
nbs1=find(squeeze(ov>quantile(ov,.9)));
%quantile(ov,.99)
%sum(visited(nbs)==0)
%sum(visited(nbs1)>0)/length(nbs1)
%parent
%y=zeros(1,n);
if length(y)==0
    y=Y(parent,:);
end
for j=nbs1'
    % visited(j)
    if visited(j)>0||visited(j)<-5
        continue
    end
    
    
    %overlap=find(N(j,:)>0&N(parent,:)>0);
    %keyboard
    overlap=find(N(j,:)>0&sum(N(visited==1,:),1)>0);
    
    S= sum(N(visited==1,:),1)>0;
    
    %[~,a]=max(histc(full(Y(:,S)),1:k));
    if sum(visited==1)>1
        [~,a]=max(histc(full(Y(visited==1,:)),1:k));
    else
        a=Y(visited==1,:);
    end
    %sum(visited==1)
    %size(full(Y(visited==1,:)))
    %histc(full(Y(visited==1,:)),1:k)
    a1=zeros(1,n);
    %[length(a) length(a1) max(find(S==1))]
    a1(S==1)=a(S==1);
    
    if length(unique(Y(j,overlap)))~=length(unique(a1(overlap)))
        histc(full(Y(j,overlap)),1:k)
        histc(a1(overlap),1:k)
        %   [length(unique(Y(j,overlap))) length(unique(a1(overlap)))]
        disp('fewer labels in one cluster')
    end
    conf=myconfusionmat(Y(j,overlap),a1(overlap),k,n);
    %keyboard
    
    
    
    %conf1
    % conf=myconfusionmat(Y(j,:),a1,k,n);
    %conf
    %pause
    % [perm]=bipartite_matching_exhaustive(conf);
    [perm]=bipartite_matching2(conf);
    
    
    % perm
    Y(j,N(j,:)>0)=perm(Y(j,N(j,:)>0));
    %obj=full(nmi(y(overlap),Y(parent,overlap)));
    %myconfusionmat(Y(j,overlap),a1(overlap),k,n);
    obj=sum(a1(overlap)==Y(j,overlap))/length(overlap);
    %perm'
    %keyboard
    %obj=1-sum(diag(conf))/sum(conf(:));
    score(j)=obj;
    %keyboard
    
    if obj<thresh || ~isreal(obj) || isnan(obj)
        %visited(j)=visited(j)-1;
        visited(j)=2;
        if verbose
            fprintf(1,'\n Horrid score %f; moving on', full(obj));
        end
        continue
    end
    
    %disp('Ha')
    
    if (0&&verbose)
        j
        conf
        obj
        %plot(a1,'.')
        imagesc(Y)
        pause
        
        %keyboard
    end
    if 0&&Y(j,10)>0
        perm'
        fprintf(1,'%d ',full(Y(j,10)));
    end
    
    [Y,N,parent,visited,score,y]=consistent_sequentially_patching(N,Y,j,visited,score,k,y);
    
    
    y=zeros(1,n);
    %y=zeros(n,k);
    a=sum(N(visited==1,:));
    tau=quantile(a,.01);
    
    y(a>tau)=myround(full(Y),visited,a,tau,k);
    
    %y(a>tau)=round(sum(Y(visited==1,a>tau))./a(a>tau));
    %y(a>tau)=round(sum(Y(visited==1,a>tau))./sum(N(visited==1,a>tau)));
    %y=kmeans(sum(Y(visited==1,:))./sum(N(visited==1,:)),k);
    
end

%end
%plot(sum(Y(visited==1,:))./sum(N(visited==1,:)));
%keyboard
end

