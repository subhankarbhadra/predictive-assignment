function [score,all_hist,visited]=consistent_sequentially_patching_faster1(N,Y,path,k,verbose)
thresh=.55;
num_seeds=size(N,1);
visited=zeros(1,num_seeds);
score=visited;

n=size(Y,2);
all_hist=zeros(k,n);
if nargin==7    
    verbose=0;
end
%first=N{1};t=sparse(n,1);t(first)=1;


parent=path(1);

%for i=1:num_seeds
%visited(parent)
if visited(parent)>0
    %disp('blah')
    return;
end

visited(parent)=1;
%[~,inds]=sort(N*N(parent,:)','descend');
%nbs1=inds;
%ov=N*N(parent,:)';
%median(ov)
%nbs1=find(squeeze(ov>quantile(ov,.9)));

%quantile(ov,.99)
%sum(visited(nbs)==0)
%sum(visited(nbs1)>0)/length(nbs1)
%parent
%y=zeros(1,n);

visited_stat = sum(N(visited==1,:),1);
S= visited_stat>0;
Sinds=find(S);

%y=Y(parent,:);
new_labels=Y(parent,Sinds);
indexes=sub2ind(size(all_hist),new_labels',Sinds');
all_hist(indexes)=1;


for j=path(2:end)
    if visited(j)>0
        continue;
    else
        visited(j)=1;
    end
    this_nbhood=N(j,:)>0;
    overlap=find(this_nbhood&S==1);
     if sum(visited==1)>1
        %[~,a]=max(histc(full(Y(visited==1,:)),1:k));
        %a=mode(Y(visited==1,:));
        
        [~,a]=max(all_hist);
     %   keyboard
    else
        a=Y(visited==1,:);
    end
    %sum(visited==1)
    %size(full(Y(visited==1,:)))
    %histc(full(Y(visited==1,:)),1:k)
    a1=zeros(1,n);
    %[length(a) length(a1) max(find(S==1))]
    %a1(S==1)=a(S==1);
    a1(Sinds)=a(Sinds);
    if 0&&length(unique(Y(j,overlap)))~=length(unique(a1(overlap)))
        histc(full(Y(j,overlap)),1:k)
        histc(a1(overlap),1:k)
        %   [length(unique(Y(j,overlap))) length(unique(a1(overlap)))]
        disp('fewer labels in one cluster')
    end
    conf=myconfusionmat(Y(j,overlap),a1(overlap),k,n);
    
    [perm]=bipartite_matching2(conf);
    
    
    % perm
    this_nbhood_indset=find(this_nbhood);
    new_labels=perm(Y(j,this_nbhood_indset));
    Y(j,this_nbhood_indset)=new_labels;
    
    indexes=sub2ind(size(all_hist),new_labels',this_nbhood_indset);
    all_hist(indexes)=all_hist(indexes)+1;
  %  Sinds
   % conf
   %keyboard
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
    
   % [Y,N,parent,visited,score,y,all_hist]=consistent_sequentially_patching_faster(N,Y,j,visited,score,k,y,0,adj_graph,all_hist);
    
    
    %y=sparse(zeros(1,n));
end


end

