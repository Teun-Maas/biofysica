function p = gumbelfun(x,mu,omega,alpha)
% P = GUMBELFUN(X,MU,OMEGA,ALPHA)
%
% The Gumbel function, parametrized by
% mean MU, width OMEGA, and probabilities that determine the width ALPHA
% and 1-ALPHA
if nargin<4
	alpha = 0.1;
end
z1		= log(-log(alpha));
z2		= log(-log(1-alpha));
z		= log(-log(0.5));
p		= 1-exp( -exp( (z1-z2)/omega * (x-mu) + z ) );
