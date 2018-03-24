function Corr=radj(corr,N,p)

% R = radj(R,N,p)
%
% Adjust Correlation Coefficient, R, for number of parameters
% N = sample size
% p = number of parameters
%
% MW

Corr= sqrt( 1-(1-corr^2)*(N-1)/(N-p-1) );