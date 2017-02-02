function p = weibullfun(x,mu,slope,~)
% P = WEIBULLFUN(X,MU,OMEGA,ALPHA)
%
% The Weibull function, parametrized by
% mean MU, slope OMEGA
p = 1-exp( -exp( 2*slope*mu/log(2) * (log(x)-log(mu))+log(log(2)) ) );
