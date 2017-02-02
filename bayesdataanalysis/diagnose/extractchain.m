function samples = extractchain(samples)
% MCMC = EXTRACTCHAIN(MCMC)
%
% change the fields in the MCMC structure to contain only 1- or 2-D
% matrices
%
% See also DIAGMCMC

[nPar,parNames] = nvar(samples); % get all parameter names
for ii			= 1:nPar
	n			= size(samples.(parNames{ii}));
	switch numel(n)
		case 2
			samples.(parNames{ii}) = samples.(parNames{ii})(:);
		case 3
			samples.(parNames{ii}) = reshape(samples.(parNames{ii}),n(1)*n(2),n(3));
		case 4
			samples.(parNames{ii}) = reshape(samples.(parNames{ii}),n(1)*n(2),n(3),n(4));
	end
end