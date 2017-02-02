function chi = probit(cdf)
% CHI = PROBIT(CDF)
%
% The probit function is the quantile function, i.e., the inverse
% cumulative distribution function (CDF), associated with the standard
% normal distribution. 
%
% This is useful for plotting reaction times.
%

myerf       = 2*cdf - 1;
myerfinv    = sqrt(2)*erfinv(myerf);
chi         = myerfinv;