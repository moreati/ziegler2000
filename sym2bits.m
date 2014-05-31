function [bits] = sym2bits(symbols)

[r,c] = size(symbols)
if c>1
   b = de2bi(symbols);
   [r2,c2] = size(b);
   bits = [];
   for i=1:c
      x = de2bi(symbols(:,i),c2);
      bits = [bits x];
   end
else
   bits = de2bi(symbols);
end




