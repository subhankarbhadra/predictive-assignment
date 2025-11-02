function C=myconfusionmat(class1,class2,k,n)
Z1=zeros(n,k);
Z2=zeros(n,k);

for i=1:k
    Z1(class1==i,i)=1;
    Z2(class2==i,i)=1;
end
C=Z1'*Z2;
if(0)
    C=zeros(k,k);
    
    
    C1=confusionmat(class1,class2);
    if size(C1,1)<k || size(C1,2)<k
        a1=unique(class1);a2=unique(class2);
        a1
        a2
        C1
        C(a1,a2)=C1;
    else
        C=confusionmat(class1,class2);
    end
    if size(C,1)<k
        keyboard
    end
end
end