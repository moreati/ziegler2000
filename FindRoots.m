function [ExtFieldEst,roots,CWLenEst] = ...
         FindRoots(FileName,WordCount,q,n)

% Look for roots in extension fields assuming base field is nonbinary

BitsPerSym = log2(q);
nmin = 6;
nmax = 255;

EncodedBitCount = nmax*WordCount*BitsPerSym;
s = sprintf('%s.fec',FileName);
CodedFile = fopen(s,'r');
[r,WordsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
if(WordsRead ~= EncodedBitCount)
   disp 'Error reading file'
   WordsRead
end

% Search fields covering all subfields from GF(2) to GF(2^20)
m = [12 20 14 16 18 11 15];
n2 = [];
ExtFields = [];

for mi=1:length(m)
   qm = 2.^m(mi)
   % Form array of valid block lengths for this field
   if (nargin < 4),
      n = SearchLens(qm-1,nmin,nmax)
   end
   n1 = [];
   MaxRoots = 0;
   for ni=1:length(n)
      fprintf('Searching for roots in GF(%d), n = %d\n',qm,n(ni));
      if BitsPerSym > 1
         EncodedBitCount = n(ni)*WordCount*BitsPerSym;
         r2 = reshape(r(1:EncodedBitCount),BitsPerSym,...
                       EncodedBitCount/BitsPerSym)';
        r2 = bits2sym(r2,BitsPerSym,EncodedBitCount/BitsPerSym);
        r2 = reshape(r2,n(ni),WordCount)';
     else
        EncodedBitCount = n(ni)*WordCount;
        r2 = reshape(r(1:EncodedBitCount),n(ni),WordCount)';
     end

     % Check combined spectrum of 3 code words for zeros
     R = CumGFFT2(r2,3,n(ni),q,qm,1);
     roots = find(~R);    % find spectral zeros
     if length(roots) > 0
        % This length has passed the 1st root test
        n1 = [n1 n(ni)];
        % Check for spectral zeros over many code words
        R = CumGFFT2(r2,WordCount,n(ni),q,qm,1);
        roots = find(~R)
        if length(roots) > 0
           % This length has passed the 2nd root test
           n2 = [n2 n(ni)]
           ExtFields = [ExtFields qm];
           break
        end
     end
  end      % for ni=1:length(n)

   if length(n2) > 0
      % there was a winner in this field, so skip the other fields
      break
   end
end         % for mi=1:length(m)


if length(n2) > 0
   % Find smallest base field with roots, and search any subfields
   % until the smallest field containing roots is found
   [EF,EFI] = sort(ExtFields);
   qm = EF(1);

  switch qm               %   search subfields of q^m
     case 4096,           %   m = 12 = 2*2*3
        m = [3 4 6 12];   %   note: GF(2) and GF(4) have no valid lengths
     case 1048576,        %   m = 20 = 2*2*5
        m = [5 10 20];
     case 16384,          % m = 14 = 2*7
        m = [7 14];
     case 65536,          % m = 16 = 2*2*2*2
        m = [8 16];
     case 262144,         % m = 18 = 2*3*3
        m = [9 18];
     case 2048,
        m = 11;
     case 32768,
        m = 15;
     otherwise,
        %return
  end

else
   m = [3:12,14:16,18,20];    % search all fields to cover RS codes
   % IMPORTANT NOTE: FIELDS GF(q^m) FOR WHICH m = 13, 17, AND 19
   % HAVE NO VALID LENGTHS n, AND THEREFORE ARE NOT REQUIRED!
   % (e.g., 8K, 128K, 512K)
end

n2 = [];


for mi=1:length(m)
   qm = 2.^m(mi);
   if rem(qm-1,q-1)<1            % bypass invalid nonbinary subfields
      % Form array of valid block lengths for this field
      if (nargin < 4),
         n = SearchLens(qm-1,nmin,nmax)
      end
      n1 = [];
      MaxRoots = 0;
      for ni=1:length(n)
         fprintf('Searching for roots in GF(%d), n = %d\n',qm,n(ni));
         if BitsPerSym > 1
            EncodedBitCount = n(ni)*WordCount*BitsPerSym;
            r2 = reshape(r(1:EncodedBitCount),BitsPerSym,...
                           EncodedBitCount/BitsPerSym)';
            r2 = bits2sym(r2,BitsPerSym,EncodedBitCount/BitsPerSym);
            r2 = reshape(r2,n(ni),WordCount)';
         else
            EncodedBitCount = n(ni)*WordCount;
            r2 = reshape(r(1:EncodedBitCount),n(ni),WordCount)';
         end
         % Check combined spectrum of 3 code words for zeros
         R = CumGFFT2(r2,3,n(ni),q,qm,1);
         roots = find(~R);    % find spectral zeros
         if length(roots) > 0
            % This length has passed the 1st root test
            n1 = [n1 n(ni)];
            % Check for spectral zeros over many code words
            R = CumGFFT2(r2,WordCount,n(ni),q,qm,1)
            roots = find(~R)
            if length(roots) > 0
               % This length has passed the 2nd root test
               n2 = n(ni)
               ExtFields = qm;
               break
            end
         end
      end      % for ni=1:length(n)

     if length(n2) > 0
        % There was a winner in this field, so skip the other fields
        break
     end
  end      % if rem(qm1,q-1)<1
end        % for mi=1:length(m)

fclose(CodedFile);
CWLenEst = n2;
ExtFieldEst = ExtFields;

if length(n2) < 1
   disp 'NOT A CYCLIC CODE'
end

