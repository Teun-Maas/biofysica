function samples = stan2jagsmat(fit)
% SAMPLES = STAN2JAGSMAT(FIT)
%
% Extract sample structure similar to MATJAGS from STAN fit object
%
% See also STAN, MATJAGS

samplesstruct = fit.extract('permuted',false,'par',{'theta','omega','lambda','gamma','g'});

parameterNames		= fieldnames(samplesstruct); % get all parameter names
for ii				= 1:numel(parameterNames)
	samples.(parameterNames{ii}) = [samplesstruct.(parameterNames{ii})]';
end

if isfield(samples,'g')
	samples.gamma = samples.g;
	samples = rmfield(samples,'g');
end

