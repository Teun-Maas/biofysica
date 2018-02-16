function phi = probit(x)
% PHI = PROBIT(X)
%
% The probit function is the quantile function, i.e., the inverse
% cumulative distribution function (X), associated with the standard
% normal distribution. 
%
% This is useful for plotting reaction times.
%

phi    = sqrt(2)*erfinv(2*x - 1);


