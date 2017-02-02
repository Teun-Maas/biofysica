function p = probitfun(x,mu,omega,alpha)
% P = PROBITFUN(X,MU,OMEGA,ALPHA)
%
% The cdf of the normal distribution, parametrized by
% mean MU, width OMEGA, and probabilities that determine the width ALPHA
% and 1-ALPHA
if nargin<4
	alpha = 0.1;
end
z		= norminv(1-alpha) -norminv(alpha);
sigma	= omega/z;
p		= normcdf(x,mu,sigma);
