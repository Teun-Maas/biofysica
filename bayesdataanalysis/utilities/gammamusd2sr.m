function [shape,rate] = gammamusd2sr(mu,sd)
% [SHAPE,RATE] = GAMMAMUSD2SR(MU,SD)
%
% Get SHAPE and RATE parameters of gamma function from mean MU and standard deviation SD
%
% SHAPE = MU^2/SD^2
% RATE = MU/SD^2
%
% Note that RATE = 1/SCALE (for Matlab's gamma functions)
%
% See also GAMMASR2MUSD, GAMFIT, GAMPDF

if mu<=0
	error('MU must be larger than 0.')
end
if sd<=0
	error('SD must be >larger than 0.')
end

shape	= mu.^2./sd.^2;
rate	= mu./sd.^2;

