function [y] = power2symGFx(x)

% Convert powers {-inf=0,0=1,1=alpha,2=alpha^2,...}
% to symbols {0,1,2=alpha,3=alpha^2,...}
y = log2(2.^(x+1) + (x<0));




