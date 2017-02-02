function samples = thinchain(samples,thinStep)
% MCMC = THINCHAIN(MCMC)
%
% thin the chains in the MCMC structure 
%
% See also DIAGMCMC, EXTRACTCHAIN

[nPar,parNames] = nvar(samples); % get all parameter names
for ii			= 1:nPar
	M			= samples.(parNames{ii}); % Matrix
	n			= size(M);
	switch numel(n)
		case 2
			samples.(parNames{ii}) = M(:,1:thinStep:end);
		case 3
			samples.(parNames{ii}) = M(:,1:thinStep:end,:);
	end
end