function [mu,sd] = gammasr2musd(shape,rate)
% [MU,SD] = GAMMASR2MUSD(SHAPE,RATE)
%
% Get mean and standard deviation parameters of gamma function from SHAPE and RATE
%
% MU = SHAPE/RATE
% SD = SQRT(SHAPE)/RATE
%
% Note that RATE = 1/SCALE (for Matlab's gamma functions)
%
% See also GAMMAMUSD2SR, GAMFIT, GAMPDF

mu = shape./rate;
sd = sqrt(shape)./rate;
