function p = logisticfun(x,theta,omega,alpha)
% P = LOGISTICFUN(X,THETA,OMEGA,ALPHA)
%
% The logistic function (logit response function in GLMs), parametrized by
% threshold THETA, width OMEGA, and probabilities that determine the width ALPHA
% and 1-ALPHA
%
% See also PSIFIT, PSIFUN, PROBITFUN, GUMBELFUN, REVGUMBELFUN,
% WEIBULLFUN, REVWEIBULLFUN
if nargin<4
	alpha = 0.1;
end
z = 2*log(1/alpha-1);
p = 1./(1+exp(-z./omega*(x-theta)));