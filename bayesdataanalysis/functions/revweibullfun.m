function p = revweibullfun(x,mu,slope,~)
% P = REVWEIBULLFUN(X,MU,OMEGA,ALPHA)
%
% The reverseWeibull function, parametrized by
% mean MU, slope OMEGA
p = exp( -exp( -2*slope*mu/log(2) * (log(x)-log(mu))+log(log(2)) ) );
