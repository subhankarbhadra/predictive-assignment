function comps=connectedComps(A)
n=size(A,1);
U=ones(1,n);
i=1;
while 1 
   rem=find(U);
   
   root=rem(1);
   v=sparse(1,n);v(root)=1;
   GO=1;
   u=v;
   old = sparse(1,n);
   while GO
       v=v*A;
       u=u+v;
       u(u>0)=1;
       v(v>0)=1;
       if (full(sum(u))==full(sum(old)))
           break;
       end
       old=u;
   end
   U(u>0)=0;
   comps{i}=find(u);
   i=i+1;
   if isempty(find(U))
       break;
   end
end  

end