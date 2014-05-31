function [z] = SumGFx(x)

% Sum the elements of vector x using GF(q) vector addition
z = x(1);
for i=2:length(x)
   z = bitxor(z,x(i));
end




