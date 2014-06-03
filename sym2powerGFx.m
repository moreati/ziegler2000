function [y] = sym2powerGFx(x)

% Convert symbols {0,1,2=alpha,3=alpha^2,...}
% to powers {-inf=0,0=1,1=alpha,2=alpha^2,...}

y = log2(~~x) + x-1;

