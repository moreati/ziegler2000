function [ConvCWLenEst,InputLenEst,ConstraintLenEst] = ...
         FindMem6(FileName,kmax,nmax,Kmax)

% RREF TECHNIQUE IN GF(2) ADAPTED TO CONVOLUTIONAL CODES

% Increase max memory length to ensure valid test up to Kmax;
% this avoids false detection of block codes
margin = 4;
Kmax = Kmax+margin;
InputLenEst = 0;
ConvCWLenEst = 0;
ConstraintLenEst = 0;
compr = ones(nmax,Kmax);
ConsecCompr = 0;

for n=2:nmax
   n
   k = 1;
   TH = (k+1)/n;
   ConstraintLenEst = 0;
   PrevRate = 1;
   PrevK = 2;
   PrevDelta = 0;
   FirstCompr = 0;

  for K=3:Kmax
     nK = n*K;
     WordCount = 2*nK;
     EncodedBitCount = nK*WordCount;
     s = sprintf('%s.fec',FileName);
     CodedFile = fopen(s,'r');

     [r,WordsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
     if(WordsRead ~= EncodedBitCount)
        disp 'Error reading file'
        WordsRead
     end

     r2 = reshape(r,nK,WordCount)';
     %[G2,nK] = RrefGF2(r2);
     [G2,kK] = RrefGF2c(r2);    % MEX C-code implementation

     rate = kK/nK;
     compr(n,K) = rate;   % DEBUG: bitwise RREF compression
     delta = nK-kK;
     delta2 = delta - PrevDelta;
     PrevDelta = delta;


     if FirstCompr
        k = n-delta2;
        TH = (k+1)/n;
        if TH > 1
           TH = 1;
        end
     end

     if rate < 1
        fprintf('Rate = %d/%d = %f, n = %d, K = %d, delta = %d, ...
                 delta2 = %d\n',kK,nK,rate,n,K,delta,delta2)
        FirstCompr = 1;

         if (rate < PrevRate) & (K == PrevK+1)
            PrevRate = rate;
            PrevK = K;
            ConvCWLenEst = n;
            ConsecCompr = ConsecCompr+1;
            if (rate < TH) & (ConstraintLenEst < 1)
               InputLenEst = k
               % Set K at first rate to drop below the threshold
               ConstraintLenEst = ((K-1)/k)+1;
               fprintf('(Rate = %f) < (threshold = %f) at ...
                        K = %d\n',rate,TH,K)
            end
         else
            % Invalidate this trial since its rate is not
            % monotonically decreasing
            ConvCWLenEst = 0;
            InputLenEst = 0;
            break
         end
      elseif ~FirstCompr
         PrevK = K;
      else
         % Invalidate this trial since its rate is not
         % monotonically decreasing
         ConvCWLenEst = 0;
         InputLenEst = 0;
         break
      end

     fclose(CodedFile);
  end

  if (ConvCWLenEst > 0) & (InputLenEst > 0) & (ConsecCompr >= margin)
     break
  end

end

% Repeat search for K if k > 1, given (n,k) and appropriate threshold,
% to revise estimate of constraint length
if (ConvCWLenEst > 0) & (InputLenEst > 1)
   TH = (k+1)/n;
   for K=3:Kmax
      nK = n*K;
      WordCount = 2*nK;
      EncodedBitCount = nK*WordCount;
      s = sprintf('%s.fec',FileName);
      CodedFile = fopen(s,'r');
      [r,WordsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
      if(WordsRead ~= EncodedBitCount)
         disp 'Error reading file'
         WordsRead
      end
      r2 = reshape(r,nK,WordCount)';
      [G2,kK] = RrefGF2c(r2);    % MEX C-code implementation
      rate = kK/nK;

      if (rate < TH)
         % Set K at first rate to drop below the threshold
         ConstraintLenEst = ((K-1)/k)+1;
         fprintf('(Rate = %f) < (threshold = %f) at K2 = %d; ...
                  K = %d\n',rate,TH,K,ConstraintLenEst)
         fclose(CodedFile);
         break
      end
      fclose(CodedFile);
   end
end

if (ConvCWLenEst > 0) & (ConstraintLenEst > 0)
   fprintf('Found convolutional code: n = %d, k = %d, K = %d\n', ...
            ConvCWLenEst,InputLenEst,ConstraintLenEst)
else
   disp 'CONVOLUTIONAL CODE NOT FOUND'
end

if 0
figure
Ki = [3:Kmax];
plot(Ki,compr(2,Ki),'b*-')
hold
plot(Ki,compr(3,Ki),'ro-')
plot(Ki,compr(4,Ki),'gd-')
xlabel('Constraint Length (K)')
ylabel('RREF Compression')
legend('Rate 1/2','Rate 1/3','Rate 1/4',0)
end

