function samples = stan2mat(fit)
% SAMPLES = STAN2JAGSMAT(FIT)
%
% Extract sample structure similar to MATJAGS from STAN fit object
%
% See also STAN, MATJAGS

samplesstruct = fit.extract('permuted',false,'par',{'theta','omega','lambda','gamma'});

parameterNames		= fieldnames(samplesstruct); % get all parameter names
for ii				= 1:numel(parameterNames)
	samples.(parameterNames{ii}) = [samplesstruct.(parameterNames{ii})]';
end

