function [RateEst,CWLenEst,InputLenEst,AlphabetSizeEst] = ...
         EstRate3(FileName,nmin,nmax,lmax)

disp 'Searching for code rate via RREF in GF(2)'
RateEst = 1;
nqmin = nmin;
nqmax = nmax;

for n = nmin:nmax
   WordCount = 2*n;     % should be > n symbols for accurate results
   EncodedBitCount = n*WordCount;
   s = sprintf('%s.fec',FileName);
   CodedFile = fopen(s,'r');
   [r,BitsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
   fclose(CodedFile);
   if(BitsRead ~= EncodedBitCount)
      disp 'Error reading file'
      BitsRead
      return
   end

   r2 = reshape(r,n,WordCount)';
   %[G2,k] = RrefGF2(r2);
   [G2,k] = RrefGF2c(r2);     % MEX C-code implementation
   rate = k/n;

   if (rate<RateEst)
      RateEst = rate;
      CWLenEst = n;
      InputLenEst = k;
      AlphabetSizeEst = 2;
      if(rate<1)
         fprintf('Code rate estimate: %d/%d = %f\n',k,n,RateEst)
      end
   end
end


for BitsPerSym = 2:lmax
   q = 2^BitsPerSym;
   if RateEst<1
      nmin = CWLenEst*log2(AlphabetSizeEst)/BitsPerSym;
      if (abs(nmin - round(nmin)) > 1e-10) | (nmin<7)
         nmin = [];
      else
         fprintf('\nSearching for code rate via RREF in GF(%d)\n',q)
      end
      nmax = nmin;
   else
      fprintf('\nSearching for code rate via RREF in GF(%d)\n',q)
   end

  for n = nmin:nmax
     nq = BitsPerSym*n;
     WordCount = n+10;    % should be > n symbols for accurate results
     EncodedBitCount = nq*WordCount;
     s = sprintf('%s.fec',FileName);
     CodedFile = fopen(s,'r');
     [r,BitsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
     fclose(CodedFile);
     if(BitsRead ~= EncodedBitCount)
        disp 'Error reading file'
        BitsRead
        EncodedBitCount
        return
     end

     EncodedBitCount = EncodedBitCount/BitsPerSym;
     r = reshape(r,BitsPerSym,EncodedBitCount)';
     r = bits2sym(r,BitsPerSym,EncodedBitCount);
     r2 = reshape(r,n,WordCount)';
     [G2,k] = RrefGFxc(r2,q);      % MEX C-code implementation
     rate = k/n;

      if (rate<=RateEst)
         RateEst = rate;
         CWLenEst = n;
         InputLenEst = k;
         AlphabetSizeEst = q;
         if(rate<1)
            fprintf('Code rate estimate: %d/%d = %f\n',k,n,RateEst)
         end
         if(rate<0.97)        % Max BCH rate = 247/255 = .969
            break             % STOP AT FIRST CLEAR SIGN OF CODING?
         end
      end
   end
end

if RateEst==1
   disp 'INPUT IS NOT CODED'
   return
end

