function [m,confM,label,label2]= cluster_acc(idx,trueIdx)
%function m= cluster_acc(idx,trueIdx)

if max(trueIdx)>5
   %disp('K too big to do all permutations. Returning 1-nmi instead'); 
   confM=confusionmat(trueIdx,idx);
   m=1-nmi(idx,trueIdx);
   return
end

v=unique(idx);
v1=unique(trueIdx);

%check if any class label is completely missing from the trueIdx
tmp=trueIdx;
if max(trueIdx)~=length(v1)
    
    for i=1:length(v1)
        tmp(trueIdx==v1(i))=i;
    end
    %[length(tmp) length(idx)]
    trueIdx=tmp;
    
end
%check if any class label is completely missing from the estimated idx
tmp=idx;
if max(idx)~=length(v)
    
    for i=1:length(v)
        tmp(idx==v(i))=i;
    end
    %[length(tmp) length(idx)]
   idx=tmp;
    
end
v=unique(idx);
v1=unique(trueIdx);

if length(v)==2&&length(v1)==2
%if length(v1)==2
    m=cluster_acc_twoclass(idx,trueIdx);
    confM=confusionmat(trueIdx,idx);
    
else
    u=perms(1:max(max(v1),max(v)));
    idxM=zeros(length(u),length(idx));
    for i=1:size(u,1)
        idx1=zeros(size(trueIdx));
        for j=1:length(v)
            idx1(idx==j)=u(i,j);
            
        end
        idxM(i,:)=idx1;
        err(i)=sum(trueIdx~=idx1);
        
    end
    
    [m,indm]=min(err);
    
    %error('This is only for two clusters');
    m=m/length(idx);
    confM=confusionmat(trueIdx,idxM(indm,:));
    label=trueIdx~=idxM(indm,:)';
    label2=idxM(indm,:);
    %end
    %confM
    if m>.6
        %keyboard;
    end
    %keyboard
end
end


function acc=cluster_acc_twoclass(idx,trueIdx)
m=0;n=length(idx);
%keyboard
m=min(sum(idx~=trueIdx),sum(3-idx~=trueIdx));
acc=m/n;

end