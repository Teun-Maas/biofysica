function BF = bayesfactor(samplesPost,samplesPrior,crit)
% Quick and dirty code
%
% BF = BAYESFACTOR(SAMPLESPOST,SAMPLESPRIOR)
%
% Determine Bayes factor for prior and posterior MCMC samples via the
% Savage-Dickey method.
%
% Determine if posterior is significantly different from a critical value
% CRIT:
% BF = BAYESFACTOR(SAMPLESPOST,SAMPLESPRIOR,CRIT)
% By default this will the null hypothesis.
%
% See also KSDENSITY

% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com


eps		= 0.01; % bin size
binse	= -100:eps:100; % bin range
if nargin<3
crit	= 0; % test null-hypothesis
end

% Posterior
[f,xi]		= ksdensity(samplesPost,'kernel','normal');
[~,indk]	= min(abs(xi-crit));

% Prior on Delta
tmp			= samplesPrior;
tmp			= tmp(tmp>binse(1)&tmp<binse(end));
[f2,x2]		= ksdensity(tmp,'kernel','normal','support',[binse(1) binse(end)]);
[~,indk2]	= min(abs(x2-crit));

v1 = f(indk);
v2 = f2(indk2);
BF = v1/v2;









