function y=myround(Y,visited,a,tau,k)

M=Y(visited==1,a>tau);
%denom=a(a>tau);
M1=histc(M,1:k)./(ones(k,1)*a(a>tau));
tmp=round(M1);
%S=find(sum(tmp,1)>1);

y=sum(diag(1:k)*tmp);

if sum(y>k)>0
   y(y>k)=randsample(k,sum(y>k),'true'); 
end

end