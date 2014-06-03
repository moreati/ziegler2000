function [Xcum] = CumGFFT2(x,WordCount,n,q,qm,mx)

% Cumulative Galois Field Fourier Transform:
% "OR" all frequency words together to isolate spectral zeros (roots).

X = GFFT2(x,WordCount,n,q,qm,mx);
Xcum = zeros(1,n);

for i=1:WordCount
   Xcum = bitor(Xcum,X(i,:));
end

