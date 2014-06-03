function [n,k,M,G,H] = Hamming(m)

% Hamming code w/ parameters m>=2, t=1
%   Ex: m=3 -> (7,4) code
n   = 2^m-1;
k   = n-m;
M   = 2^k;       % code dimension (# of code words)

% Generator polynomials for systematic code
P = zeros(k,n-k);
Pm = zeros(k,1);
j = 3;
ii = 2;
for i=1:k
   % Bypass all power-of-2 elements, since these
   % represent the systematic part
   if abs(j-2^ii) < 1e-10
      j = j+1;
      ii = ii+1;
   end
   Pm(i) = j;
   j = j+1;
   for l=1:(n-k)
      P(i,l) = rem(floor(Pm(i)/(2^(n-k-l))),2);
   end
end
G = [P eye(k)];            % generator matrix
H = [eye(n-k) P'];         % parity check matrix

