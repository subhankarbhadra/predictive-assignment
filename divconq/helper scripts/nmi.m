function x = nmi(C,Chat)
N = length(C);
CF = confusionmat(C,Chat);
p1 = sum(CF,1)/N;
p2 = sum(CF,2)/N;
p12 = reshape(CF,[],1)/N;
H1 = NtroP(p1);
H2 = NtroP(p2);
H12 = NtroP(p12);
MI = H1 + H2 - H12;
%x = 2*(1 - H12/(H1 + H2)); % Proposed by Danon et al. (2005), 'compare' function in the 'igraph' package in R
%uses this
x = MI/(min(H1, H2)+1e-100); % This is the definition Newman uses. 
%keyboard
end