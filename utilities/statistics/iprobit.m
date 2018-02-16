function x = iprobit(phi)
% PHI = PROBIT(X)
%
% The probit function is the quantile function, i.e., the inverse
% cumulative distribution function (X), associated with the standard
% normal distribution. 
%
% This is useful for plotting reaction times.
%


x = (erf(phi/sqrt(2))+1)/2;
