function [Gtrain,Gtest,Gtest2,idx2]=process_real_graph(Gtrain,Gtest,Gtest2,thresh)
if nargin==3
    thresh=1;
end
d=sum(Gtrain,1);
n=length(d);
idx=1:n;

idx=find(d>=thresh);

%NIPS 2015 uses thresh=1
%idx=find(d>0); %This is used for annals for all real graphs.
%idx = find(d>1); %new annals results (it was 2, but that must be set
%after I got the results)
%idx = find(d>2);


comps=connectedComps(Gtrain(idx,idx));
for i=1:length(comps)
    s(i)=length(comps{i});
end
[~,largest]=max(s);
%largest
idx2=idx(comps{largest});
Gtrain=Gtrain(idx2,idx2);
Gtest=Gtest(idx2,idx2);
Gtest2=Gtest2(idx2,idx2);

end
