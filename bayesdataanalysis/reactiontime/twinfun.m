function RT = twinfun(T,lambda1,lambda2,mu,omega,delta,gamma,kappa)
% RT = TWINFUN(T,LAMBDA1,LAMBDA2,OMEGA,DELTA,GAMMA,KAPPA)
%
% See also TWINPI, TWINPW

% Probability of integration
PI = twinpi(T,lambda1,lambda2,omega);
% Probability of warning
PW = twinpw(T,lambda1,lambda2,gamma);
% Observed mean reaction time
RT = 1./lambda1+mu-delta*PI-kappa*PW;



