function [n,k,M,G,H] = Golay23(gen)

%   Binary (23,12,7) Golay code
%   dmin=7, t=3
n   = 23;
k   = 12;
M   = 2^k;       % code dimension (# of code words)

if gen == 1
   % Generator polynomials for general code
   g = [1 1 0 0 0 1 1 1 0 1 0 1];
else
   g = [1 0 1 0 1 1 1 0 0 0 1 1];
end

G = zeros(k,n);
for i=1:k
   G(i,:) = [zeros(1,i-1) g zeros(1,n-k-i+1)];   % generator matrix
end




