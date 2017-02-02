function nChain = nchain(mcmc)
% N = NCHAIN(X)
%
% gives number of chains of an MCMC array
%
% See also NITER, NVAR


if isstruct(mcmc)
	parNames	= fieldnames(mcmc);
	mcmc = mcmc.(parNames{1}); % all parameters should have same number of iterations and chains
end

nd = ndims(mcmc);
switch nd
	case 2
		[nIter,nChain] = size(mcmc); % first is number of iterations, 2nd = number of chains
	case 3
		[nIter,nChain,~] = size(mcmc); % 3rd is number of subjects/conditions
end

if nChain>nIter
	nChain = nIter;
end

