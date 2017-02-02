function p = normpdfun(x,mu,sigma)
% P = NORMPDFUN(X,MU,SIGMA)
%
% The Normal probability density function, parametrized by
% mean MU and standard deviation SIGMA

a = 1./(sigma*sqrt(2*pi));
b = mu;
c = sigma;

p = a*exp(-(x-b).^2/(2*c.^2));