function [m,c,WordsWritten] = ...
         BlockEnc(FileName,WordCount,BlockCount,n,k,q,G,FileIO)

% Generate Galois Field vector representations for GF(2^3) to GF(2^20)
if ~exist('vectors8_m1')
   ReadFields
end
[vectors,symbols] = SelField(q,1);

WordsWritten = 0;
SourceBitCount = k*WordCount*log2(q);
EncodedBitCount = n*WordCount*log2(q);
m = zeros(WordCount,k);       % input message words
c = zeros(WordCount,n);       % output code words

if FileIO
   s = sprintf('%s.src',FileName)
   SourceFile = fopen(s,'w');
   s = sprintf('%s.fec',FileName)
   CodedFile = fopen(s,'w');
end

% Encode a random bit stream (represented as a binary matrix)
for j=1:BlockCount
   m = floor(q*rand(WordCount,k));
   if q==2
      for i=1:WordCount
         c(i,:) = rem(m(i,:)*G,q);
      end
   else
      % Convert symbols {0,1,2=alpha,...}
      % to powers {-inf=0,0=1,1=alpha,...}
      m2 = log2(~~m) + m-1;
      G2 = log2(~~G) + G-1;
      for l=1:WordCount
         for i=1:n
            a = zeros(1,k);
            for j=1:k
               ij = m2(l,j)+G2(j,i);
               if ij >= 0 % only add if input is nonzero (power > -inf)
                  ij = rem(ij,q-1);
                  a(j) = vectors(ij+2);
               else
                  a(j) = 0;
               end
            end
            cij = SumGFx(a);
            c(l,i) = symbols(cij+1);
         end
      end
   end
   m1 = fliplr(de2bi(m'));
   c1 = fliplr(de2bi(c'));
   d = reshape(m1',1,SourceBitCount);     % source bit stream
   s = reshape(c1',1,EncodedBitCount);    % Tx bit stream
   if FileIO
      fwrite(SourceFile,d,'ubit1');
      WordsWritten = WordsWritten + fwrite(CodedFile,s,'ubit1');
   end
end

if FileIO
   fclose(SourceFile);
   fclose(CodedFile);
else
   WordsWritten = 0;
end

