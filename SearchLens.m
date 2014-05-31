function [n] = SearchLens(qm1,nmin,nmax)

r   =   rem(qm1,[1:nmax]);
n   =   find(~r);
i   =   find(n>nmin);
n   =   n(i);




