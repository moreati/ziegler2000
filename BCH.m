function [n,k,M,q,G] = BCH(n,k,t,q,g)

% Check that degree of generator polynomial equals the redundancy
if length(g) ~= (n-k+1)
   disp 'Invalid BCH code'
   return
end

M = q^k;       % code dimension (# of code words)
G = zeros(k,n);
for i=1:k
   G(i,:) = [zeros(1,i-1) g zeros(1,k-i)];   % generator matrix
end

