% FEC CODE RECOGNITION MAIN PROGRAM

% Select test code to identify
%FileName = 'Uncoded'; % control test case
%%% Binary codes
%FileName = 'Hamming3';
%FileName = 'Hamming4';
%FileName = 'Golay23-1';
%FileName = 'Golay23-2';
%FileName = 'Crc7-4-1';
%FileName = 'Bch63-45-3';
%FileName = 'rm16-5-3';
%%% 4-ary codes
%FileName = 'Bch21-15-1-4';
%FileName = 'Bch21-12-2-4';
%%% 16-ary codes
%FileName = 'Bch85-73-3-16';
%%% 64-ary codes
%FileName = 'Rs63-57-3';
%%% 256-ary codes
%FileName = 'Rs255-245-5';
%%% Rate 1/n convolutional codes
%FileName = 'ConvR1-2_K3';
%FileName = 'ConvR1-2_K4';
%FileName = 'ConvR1-2_K5';
%FileName = 'ConvR1-2_K6';
%FileName = 'ConvR1-2_K7';
%FileName = 'ConvR1-2_K8';
%FileName = 'ConvR1-2_K9';
%FileName = 'ConvR1-2_K9b';
%FileName = 'ConvR1-2_K14';
%FileName = 'ConvR1-3_K3';
%FileName = 'ConvR1-3_K4b';
%FileName = 'ConvR1-3_K8';
%FileName = 'ConvR1-3_K14';
%FileName = 'ConvR1-4_K3';
%FileName = 'ConvR1-4_K7';
%FileName = 'ConvR1-4_K14';
%FileName = 'ConvR2-3_K3_M4';
%FileName = 'ConvR2-3_K6_M10';
FileName = 'ConvR3-4_K4_M9';

% Search for rate k/n convolutional codes via RREF in GF(2), including
% constraint length
[ConvCWLenEst,InputLenEst,ConstraintLenEst] = FindMem6(FileName,3,4,14);

if (ConvCWLenEst > 0) & (ConstraintLenEst > 0)
   if (InputLenEst == 1)
      % Search for generator polynomials given K for rate 1/n
      % convolutional codes
      G = FindGenPolys(FileName,ConvCWLenEst,ConstraintLenEst);
   else
      % Generator polynomial search for k>1 convolutional codes is too
      % complex
      G = [];
   end

   fprintf('\n----------------\n')
   fprintf('Identified (%d,%d,%d) convolutional code\n',...
               ConvCWLenEst,InputLenEst,ConstraintLenEst);
   if ~isempty(G)
      fprintf('Order-%d generator polynomials (in octal):\n',...
               ConstraintLenEst-1);
      G
   end
   return
end

% Generate Galois Field vector representations for GF(2^3) to GF(2^13)
if ~exist('vectors8_m1')
   ReadFields
end

WordCount = 30;
warning off         % disable log of zero warning message

% Estimate block code params (n,k,q) via RREF in GF(q)
[RateEst,CWLenEst,InputLenEst,AlphabetSizeEst] = ...
   EstRate3(FileName,7,255,8)

if RateEst<1

% Search for cyclic codes via GFFT root search in extension fields
[ExtFieldEst,ExtRoots,CWLenEst2] = ...
   FindRoots(FileName,WordCount,AlphabetSizeEst,CWLenEst)
warning on

if ~isempty(ExtFieldEst)
   if CWLenEst ~= CWLenEst2
      disp 'Error: conflicting code word lengths'
      dbcont   % exit
   end
end

if length(ExtRoots)>0
   disp 'Found cyclic block code'

  % Convert roots to base field by multiplying by (q^m-1)/n
ExtRoots = ExtRoots-1;
BaseRoots = ExtRoots*(ExtFieldEst-1)/CWLenEst

% Rebuild generator polynomial by convolving all roots together:
% g(x) = (x+alpha^r1)*(x+alpha^r2)...*(x+alpha^ri)
g = PolyGFx(BaseRoots,CWLenEst,AlphabetSizeEst,ExtFieldEst,1)

redundancy = length(g) - 1
InputLenEst2 = CWLenEst - redundancy
if InputLenEst ~= InputLenEst2
   disp 'Error: conflicting input word lengths'
   dbcont   % exit
end

[ConsecRoots,t] = FindConsecRoots(ExtRoots)
dmin = EstMinDist(FileName,CWLenEst,AlphabetSizeEst,InputLenEst*2)
if t>0
   CodeType = 'BCH'
   if (CWLenEst==AlphabetSizeEst-1 & ExtFieldEst==AlphabetSizeEst)
      CodeType = 'Reed-Solomon'
   end
   if dmin < 2*t+1
      disp 'Error: invalid minimum distance (dmin < 2*t+1)'
      dbcont   % exit
   else
      disp 'Check: dmin > 2*t'
   end
   if rem(ExtFieldEst-1,CWLenEst) == 0
      disp 'Check: q^m-1 | n'
   else
      disp 'Error: invalid BCH code (n != q^m-1)'
      dbcont   % exit
   end
else
   CodeType = 'CRC or perfect cyclic block'
   t = floor((dmin-1)/2)

  redundancy2 = EstRedundancy(CWLenEst,AlphabetSizeEst,t)
  if redundancy == redundancy2
     disp 'Assume perfect block code'
     CodeType = 'perfect block'
     if InputLenEst == 1
        if t == (CWLenEst-1)/2
           CodeType = 'binary repetition'
        else
           disp 'Error: invalid repetition code (t != (n-1)/2)'
           dbcont   % exit
        end
     elseif (InputLenEst==12 & CWLenEst==23)
        if t == 3
           CodeType = 'binary Golay'
        else
           disp 'Error: invalid Golay code (t != 3)'
           dbcont   % exit
        end
     elseif CWLenEst == ExtFieldEst-1
           if t == 1
              CodeType = 'binary Hamming'
           else
              disp 'Error: invalid Hamming code (t != 1)'
              dbcont   % exit
           end
        else
           disp 'Error: invalid perfect code (no match)'
           dbcont   % exit
        end

      else
         disp 'Assume CRC block code'
         CodeType = 'CRC block';
      end
   end
else
   disp 'Assume non-cyclic block code'
   CodeType = 'non-cyclic block';
   dmin = EstMinDist(FileName,CWLenEst,AlphabetSizeEst,InputLenEst*2)
   t = floor((dmin-1)/2)
   redundancy = CWLenEst - InputLenEst;
   redundancy2 = EstRedundancy(CWLenEst,AlphabetSizeEst,t)
   if (redundancy==redundancy2 & t==1)
      CodeType = 'binary Hamming';
   end
end

fprintf('\n----------------\n')
fprintf('Identified code:\n')
fprintf('(%d,%d) %d-ary %d-error correcting %s code\n',...
         CWLenEst,InputLenEst,AlphabetSizeEst,t,CodeType);
if ~isempty(ExtFieldEst)
   fprintf('Roots in GF(%d) as powers of alpha:\n',ExtFieldEst);
   BaseRoots
   if max(g)<2
      fprintf('Order-%d binary generator polynomial (in octal):\n', ...
               redundancy);
      g = oct2gen(g,[1,1,redundancy])
   else
      fprintf('Order-%d generator polynomial:\n',redundancy);
      g
   end
end

end

