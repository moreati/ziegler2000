function [n,k,M,q,G] = ReedMuller(r,m)

q = 2;
if r==1
   % First order Reed-Muller code w/ parameters (1,m)
   n = 2^m;
   k = m+1;
   M = 2^k;       % code dimension (# of code words)
   G = [ones(1,n); flipud(sym2bits([0:n-1]')')];   % generator matrix
else
   disp 'Higher order Reed-Muller codes not implemented'
end

