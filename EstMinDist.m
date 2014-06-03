function [dmin] = EstMinDist(FileName,n,q,WordCount)

% Estimate minimum distance between code words

% Read bit stream and segment a block into a fixed # of code words
EncodedBitCount = n*WordCount;
s = sprintf('%s.fec',FileName);
CodedFile = fopen(s,'r');
[r,BitsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
if(BitsRead ~= EncodedBitCount)
   disp 'Error reading file'
   BitsRead
end
r = reshape(r,n,WordCount)';

% Bitwise XOR and sum to find distance between each pair
d = inf*ones(WordCount);
for i=1:WordCount
   for j=(i+1):WordCount
      d(i,j) = sum(rem(r(i,:)+r(j,:),q));
      if d(i,j)<1
         d(i,j) = inf;     % ignore zero distances
      end
   end
end

dmin = min(min(d));

