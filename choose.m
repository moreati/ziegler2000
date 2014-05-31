function [nck] = choose(n,k)

% n choose k, i.e., n!/(k!(n-k)!)
nck = (prod((n-k+1):n)./prod(1:k));




