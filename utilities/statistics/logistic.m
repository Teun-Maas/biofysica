function p = logistic(x)
% P = LOGISTIC(X)
%
% The logistic or sigmoid function:
% P = 1/(1+exp(-X))
%
% See also LOGIT, PROBIT

% 2013 Marc van Wanrooij

p = 1./(1+exp(-x));
