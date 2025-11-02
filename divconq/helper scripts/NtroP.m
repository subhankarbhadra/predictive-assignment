function u = NtroP(p)
 u = sum(-(p(p>0).*(log(p(p>0)))));
end