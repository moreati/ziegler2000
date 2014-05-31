function [m,c,WordsWritten] = ...
         ConvEnc3(FileName,WordCount,n,k,K,G,FileIO)

R = k/n;                        % rate
WordsWritten = 0;




SourceBitCount = k*WordCount;
EncodedBitCount = n*(WordCount+K-1);
m = zeros(WordCount,k);       % input message words
c = zeros(n,WordCount+K-1);

if FileIO
   s = sprintf('%s.src',FileName)
   SourceFile = fopen(s,'w');
   s = sprintf('%s.fec',FileName)
   CodedFile = fopen(s,'w');
end

% Encode a random bit stream (represented as a binary matrix)
m = floor(2*rand(WordCount,k));

for i1=1:n
   for i2=1:k
      i3 = n*(i2-1)+i1;
      ci = gfconv(m(:,i2),G(i3,:));
      l = length(ci);
      % Conv function does not output entire sequence
      % if ending in a string of zeros
      c(i1,1:l) = rem(c(i1,1:l) + ci, 2);
   end
end

d = reshape(m',1,SourceBitCount);      % source bit stream
s = reshape(c,1,EncodedBitCount);      % Tx bit stream
if FileIO
   fwrite(SourceFile,d,'ubit1');
   WordsWritten = fwrite(CodedFile,s,'ubit1');
end

if FileIO
   fclose(SourceFile);
   fclose(CodedFile);
else
   WordsWritten = 0;
end
