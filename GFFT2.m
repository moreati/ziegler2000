function [X] = GFFT2(x,WordCount,n,q,qm,mx)

% Galois Field Fourier Transform:
% Time-domain vector x is over GF(q), elements are q-ary
% Frequency-domain vector X is over GF(q^m), elements are q^m-ary

[vectors,symbols] = SelField(qm,mx);
X = zeros(WordCount,n);
a = zeros(1,n);
qm1 = qm-1;
% Reverse order in time sequence to form polynomial
x = fliplr(x);

if 1
   % C code implementation
   X = gfft2c(x,vectors,symbols,q,qm);
else
   % Convert symbols {0,1,2=alpha,...} to
   % powers {-inf=0,0=1,1=alpha,...}
   x2 = sym2powerGFx(x);
   K1 = qm1/n;
   K2 = qm1/(q-1);

  for l=1:WordCount
     for j=0:(n-1)        % frequency index
        for i=0:(n-1)     % time index
           % Convert GFFT coefficients from powers of gamma to alpha:
           % Calculate gamma as an nth root of unity in GF(q^m)
           ij = rem(i*j,n);
           % Calculate alpha as an (q^m-1) root of unity in GF(q^m)
           ij = rem(ij*K1,qm1);

           % Convert input coefficients from powers of beta to alpha:
           % Calculate beta as an (q-1) root of unity in GF(q)
           k = x2(l,i+1);
           if k >= 0      % only add if input is nonzero (power > -inf)
              % Calculate alpha as an (q^m-1) root of unity in GF(q^m)
              k = rem(k*K2,qm1);
              ai = rem(ij+k,qm1);
              a(i+1) = vectors(ai+2);
           else
              a(i+1) = 0;
           end




         end
         Xij = SumGFx(a);
         X(l,j+1) = symbols(Xij+1);
      end
   end
end



