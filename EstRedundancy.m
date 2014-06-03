function redundancy = EstRedundancy(n,q,t)

Vq = 0;
for i=0:t
   Vq = Vq + choose(n,i)*(q-1)^i;
end
redundancy = log2(Vq)/log2(q);      % logq(Vq(n,t))

