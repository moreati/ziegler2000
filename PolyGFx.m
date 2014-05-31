function [v] = PolyGFx(roots,n,q,qm,mx)

% Build polynomial from Galois Field roots

[vectors,symbols] = SelField(qm,mx);
vectors = de2bi(vectors);     % gfconv requires vector array
qm1 = qm-1;
L = length(roots);
v1 = [0 roots(1)];

for l=2:L
   v2 = [0 roots(l)];
   v1 = gfconv(v1,v2,vectors);
end

% Convert coefficients back to GF(q)
v1 = v1*(q-1)/qm1;
% Convert powers {-inf=0,0=1,1=alpha,...} to symbols {0,1,2=alpha,...}
j = find(v1<0);
v = v1+1;
v(j) = 0;




