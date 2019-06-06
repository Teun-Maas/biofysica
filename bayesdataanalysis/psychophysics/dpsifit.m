function samples = dpsifit(x,y,s,e,varargin)
% MCMC = PSIFIT(X,Y,S)


%% Initialization

% check for input orientation
% Y needs to be either:
% - a nx1 column vector with 0s and 1s
% - a nx2 column matrix with rate in colum 1 and #trials in column 2
if size(y,2)>2
	warning('Y should be a 1/2-column vector/matrix');
	y = y';
	x = x';
	s = s';
end
if size(y,1)~=size(x,1)
	warning('Y and X should have same orientations');
	x = x';
end
if size(y,1)~=size(x,1)
	error('X and Y should have same number of rows');
end
if isempty(s)
	s = ones(size(y,1),1);
end
%% S should contain consecutive values, no missing values
[~,~,s] = unique(s);

% model
if size(y,2)==2
	respDist = 'binomial';
elseif size(y,2)==1
	respDist = 'bernouilli';
end


fig				= keyval('figure',varargin,200); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
alpha			= keyval('alpha',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 alpha for all subjects/conditions

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

mcmcMethod		= keyval('mcmcMethod',varargin,'jags');

% graphic flags
diagFlag		= keyval('showDiag',varargin,false); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
predFlag		= keyval('showPred',varargin,false); % show posterior predictive
centroidFlag	= keyval('showCentroid',varargin,'mean'); % mode, median, mean

%% Actual regression

	samples = jags_dpsifit(x,y,s,e,...
		'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
		'saveName',saveName,'nChains',nChains,...
		'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
		'dic',dic,'alpha',alpha,'respDist',respDist);


%% Thin
if thinSteps>1
	samples = thinchain(samples,thinSteps);
end

%% MCMC diagnostics
if diagFlag
	diagmcmc(samples)
end

%% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D







function [samples,stats] = jags_dpsifit(x,y,s,e,varargin)
% SAMPLES = JAGS_PSIFIT(X,Y)
%
% See also PSIFIT

%% Psychometric function
fun				= keyval('function',varargin,@expfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
alpha			= keyval('alpha',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 alpha for all subjects/conditions
% modelname		= fullfile(pwd, 'jags_exp_model.txt');
modelname = '/Users/marcw/Dropbox/Manuscript/Luuk van de Rijt/#4 SPIN AV CI/matlab/dpsi_model.txt';
%% MCMC parameters
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);
respDist		= keyval('respDist',varargin);
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.

%% Numbers (amount of data and groups)
Ngroups			= numel(unique(s));
Ntotal			= size(x,1);

%% Write the model
% writemodel(alpha,Ngroups,fun,respDist,modelname);

%% Data
parameters		= {'lambda','mulambda','kappalambda','omega','muomega','sigmaomega','theta','mutheta','sigmatheta','dtheta','mudtheta','sigmadtheta'};
dataStruct		= struct('x',x,'y',y,'s',s,'e',e,'Ntotal',Ntotal,'Nsubj',Ngroups); % data structure
initsStruct		= initialize(Ngroups,nChains); % initial parameter values

%% parallel?
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
% 	doparallel		= 0; % do not use parallelization

%% MCMC


fprintf( 'Running JAGS...\n' );
[samples, stats] = matjags( ...
	dataStruct, ...                     % Observed data
	modelname, ...    % File that contains model definition
	initsStruct, ...                          % Initial values for latent variables
	'doparallel' , doparallel, ...      % Parallelization flag
	'nchains', nChains,...              % Number of MCMC chains
	'nburnin', burnInSteps,...              % Number of burnin steps
	'nsamples', nIter, ...           % Number of samples to extract
	'thin', thinSteps, ...                      % Thinning parameter
	'dic',dic, ...                       % Do the DIC?
	'monitorparams', parameters, ...     % List of latent variables to monitor
	'savejagsoutput',0, ...          % Save command line output produced by JAGS?
	'verbosity',0, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
	'cleanup',1);                    % clean up of temporary files?



%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end


function initsStruct	= initialize(Nsubj,nChains)
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lapse rate) are determined from dependent variables/
% y data
for ii = 1:nChains
% 	initsStruct(ii).lambda			=zeros(Nsubj,1); % standardized mean
	initsStruct(ii).theta			=zeros(Nsubj,1); % standardized mean
end




