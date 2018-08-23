function diagmcmc(mcmcStruct,varargin)
% DIAGMCMC(MCMC)
%
% Diagnose MCMC sampling

% Utility programs for use with the book,
%    Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:
%    A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.
% This file contains several functions that are called by other programs
% or can be called directly by the user. To load all the functions into
% R's working memory, at R's command line type:

parName		=  keyval('parName',varargin);

%% If no parameter name is given, then diagnose chains from ALL parameters
if isempty(parName)
	parameterNames	= fieldnames(mcmcStruct); % get all parameter names
	for parIdx			= 1:numel(parameterNames)
		n = size(mcmcStruct.(parameterNames{parIdx}),3);
		for ii = 1:n
			figure
			a		= squeeze(mcmcStruct.(parameterNames{parIdx})(:,:,ii));
			samp	= mcmcStruct;
			samp.(parameterNames{parIdx}) = a;
			diagmcmc(samp,'parName',parameterNames{parIdx});
		end
	end
else
	%% Traces
	subplot(221)
	x = transpose(mcmcStruct.(parName));
	m = length(x);
	t = 1:m;
	plot(t,x);
	box off
	set(gca,'TickDir','out','XTick',0:1000:length(x)+100);
	xlim([-100 length(x)+100]);
	ylabel('Parameter value');
	xlabel('Iteration');
	
	%% Autocorrelation
	subplot(222)
	plotacfmcmc(mcmcStruct,'parName',parName)
	
	subplot(223)
	plotgrbmcmc(mcmcStruct,'parName',parName)
	
	subplot(224)
	plotdensmcmc(mcmcStruct','parName',parName);
end
