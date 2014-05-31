function [vectors,symbols] = GenGFx(qm,mx)

% Generate a Galois Field GF(q^m) based on a minimal polynomial m(x).
% Note that mx is provided as a binary number, not an array;
% also, mx is provided as the actual polynomial minus the x^(q^m) term,




%   since   alpha^(q^m) = mx,     m(x) = mx + x^(q^m)

vectors = zeros(qm,1);
v = 1;

for i=1:(qm-1)
   vectors(i+1) = v;
   v = v*2;
   if v>(qm-2)
      v = rem(v,qm);
      v = bitxor(v,mx);
   end
end

[y,i] = sort(vectors);
symbols = i-1;
size = sprintf('uint%d',log2(qm));
s = sprintf('gf%d_%d.vec',qm,mx);
fid = fopen(s,'w')
%fprintf(fid,'%d',vectors);
fwrite(fid,vectors,size)
fclose(fid)
s = sprintf('gf%d_%d.sym',qm,mx);
fid = fopen(s,'w')
%fprintf(fid,'%d',symbols);
fwrite(fid,symbols,size)
fclose(fid)




