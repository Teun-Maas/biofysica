% CHAIN_GLOBALS
%
% Set default number of chains and parallelness for MCMC
%
% Based on utility functions in R from:
%    Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
%    A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.


poolobj = gcp; % If no pool, do not create new one.
if isempty(poolobj)
	nCores = 1;
else
	nCores = poolobj.NumWorkers;
end

if nCores > 4,
	nChainsDefault			= 4;  % because JAGS has only 4 rng's.
	runjagsMethodDefault	= 'parallel';
elseif nCores == 4,
	nChainsDefault			= 3;  % save 1 core for other processes.
	runjagsMethodDefault	= 'parallel';
elseif nCores < 4,
	nChainsDefault			= 3;
	runjagsMethodDefault	= 'rjags'; % NOT parallel
end