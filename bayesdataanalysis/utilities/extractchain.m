function samples = extractchain(samples)
% MCMC = EXTRACTCHAIN(MCMC)
%
% Combine chain and sample columns into 1.


parameterNames	= fieldnames(samples); % get all parameter names
for ii			= 1:numel(parameterNames)
	n			= size(samples.(parameterNames{ii}));
	switch numel(n)
		case 2
			samples.(parameterNames{ii}) = samples.(parameterNames{ii})(:);
		case 3
			samples.(parameterNames{ii}) = reshape(samples.(parameterNames{ii}),n(1)*n(2),n(3));
		case 4
			samples.(parameterNames{ii}) = reshape(samples.(parameterNames{ii}),n(1)*n(2),n(3),n(4));
			
	end
end