function [vectors,symbols] = ReadGFx(qm,mx)

% Read a Galois Field GF(q^m) based on a minimal polynomial m(x).
% Note that mx is provided as a binary number, not an array;
% also, mx is provided as the actual polynomial minus the x^(q^m) term,
% since alpha^(q^m) = mx,        m(x) = mx + x^(q^m)

size = sprintf('ubit%d',log2(qm));
s = sprintf('gf%d_%d.vec',qm,mx);
fid = fopen(s,'r');
vectors = fread(fid,qm,size);
fclose(fid);

s = sprintf('gf%d_%d.sym',qm,mx);
fid = fopen(s,'r');
symbols = fread(fid,qm,size);
fclose(fid);

