function samples = regjags(y,x,varargin)
% MCMC = REGJAGS(Y,X)
%
% Simple/robust linear regression with JAGS
%
% MCMC structure contains fields with MCMC samples:
% - beta0: offset
% - beta1: slope
% - sigma: standard deviation of noise around regression line
% - nu: "normality parameter" for robust regression
% - pearson: pearson's correlation substructure.
%
% MCMC = REGJAGS(Y,X,'NAME',VALUE)
% Additional name-value pair inputs include: 
% - 'type' - 'simple' (default) or 'robust'
%
% and some JAGS sampling parameters
% - 'burnInSteps': Total number of burn-in steps
% - 'numSavedSteps': Total number of steps in chains to save
% - 'thinSteps': Number of steps to thin (to prevent autocorrelation).
% - 'saveName': save MCMCM samples in file 'saveName'
% - 'nChains': number of chains to run
% - 'runJagsMethod'
%
% You need to install JAGS and MATJAGS
%	http://mcmc-jags.sourceforge.net/
%	-  http://psiexp.ss.uci.edu/research/programs_data/jags/ and/or https://github.com/msteyvers/matjags
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij
%
% See also MATJAGS, PLOTPOST, PLOTLOC


%% Simulation
if nargin<2
	% Simulated height and weight data:
% 	[HtWtData,~,height,weight] = HtWtDataGenerator(30,5678); % from Kurschke
% 	x		= HtWtData(:,height);
% 	y		= HtWtData(:,weight);
% 	
	x = -90:90;
	x = x';
	y = 0.9*x+10*randn(size(x));
end

%% Initialization
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,500); 
thinSteps		= keyval('thinSteps',varargin,1);
burnInSteps		= keyval('burnInSteps',varargin,200);
saveName		= keyval('saveName',varargin,'HierMultLinRegressData-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
type			= keyval('type',varargin,'simple'); % simple or robust
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
diagFlag		= keyval('diag',varargin,false);
diagFlag		= keyval('diag',varargin,true);

dic				= keyval('dic',varargin,false);

%% Actual regression
[samples,stats]				= genMCMC(x,y,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,'type',type,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic);

%% diagnose
if diagFlag
	parameterNames	= fieldnames(samples); % get all parameter names
	for parIdx			= 1:numel(parameterNames)
		figure(100+parIdx)
		clf;
		diagmcmc(samples,'parName',parameterNames{parIdx});
	end
end

%%

%% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

%% Pearson correlation
pearson		= genMCMCpearson(x,y);

% samples		= addfields(samples,pearson);
samples.pearson = pearson;

%% Checks and display
posteriorprediction(x,samples);

%% Sub-functions
function [samples,stats] = genMCMC(x,y,varargin)
% B = GENMCMC(X,Y)
%
% Generate MCMC chains

numSavedSteps	= keyval('numSavedSteps',varargin);
thinSteps		= keyval('thinSteps',varargin);
burnInSteps		= keyval('burnInSteps',varargin);
saveName		= keyval('saveName',varargin);
nChains			= keyval('nChains',varargin);
type			= keyval('type',varargin);
runjagsMethod	= keyval('runjagsMethod',varargin);
dic			= keyval('dic',varargin);

% modelname = which('regress.txt');
% if ~exist(modelname,'file')
if strcmpi(type,'simple')
	writesimplemodel;
	modelname = fullfile(pwd, 'model.txt');
	parameters		= {'beta0','beta1','sigma'};		% The parameter(s) to be monitored.
elseif strcmpi(type,'robust')
	writerobustmodel;
	modelname = fullfile(pwd, 'model.txt');
	parameters		= {'beta0','beta1','sigma','nu'};		% The parameter(s) to be monitored.
	
end
% end

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
[zx,mux,sdx]		= zscore(x);
[zy,muy,sdy]		= zscore(y);
nData				= size(x,1);

% Specify data, as a structure
dataStruct = struct('x',zx,'y',zy,'Ndata',nData);

%% INTIALIZE THE CHAINS.
r			= corrcoef(x,y);
initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).beta0		= 0;    % because data are standardized
	initsStruct(ii).beta1		= r(2);        % because data are standardized
	initsStruct(ii).sigma		= (1-r(2)^2);  % because data are standardized
	if strcmpi(type','robust')
		initsStruct(ii).nu		= 100;        % because data are standardized
	end
end

%% RUN THE CHAINS
% adaptSteps		= 500;			% Number of steps to 'tune' the samplers.
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
	
else
	doparallel		= 0; % do not use parallelization
end
fprintf( 'Running JAGS...\n' );
% [samples, stats, structArray] = matjags( ...
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

%% EXAMINE THE RESULTS
% checkConvergence = false;
% if checkConvergence
% 	fnames = fieldnames(samples);
% 	n = numel(fnames)-1;
% 	for ii = 1:n
% 		% thetaSample = reshape(mcmcChain,[numSavedSteps+1 nCoins]);
% 		ns = ceil(sqrt(n));
% 		figure(1)
% 		subplot(ns,ns,ii)
% 		s		= samples.(fnames{ii});
% 		s		= s(1:16667);
% 		s		= zscore(s);
% 		maxlag	= 35;
% 		[c,l]	= xcorr(s,maxlag,'coeff');
% 		h		= stem(l(maxlag:end),c(maxlag:end),'k-');
% 		xlim([-1 maxlag]);
% 		set(h,'MarkerEdgeColor','none');
% 		axis square;
% 		box off
% 		set(gca,'TickDir','out');
% 		str = fnames{ii};
% 		xlabel('Lag');
% 		ylabel('Autocorrelation');
% 		title(str);
% 		%   show( gelman.diag( codaSamples ) )
% 		%   effectiveChainLength = effectiveSize( codaSamples )
% 		%   show( effectiveChainLength )
% 	end
% end

parameterNames	= fieldnames(samples); % get all parameter names
for ii			= 1:numel(parameterNames)
	
	z			= samples.(parameterNames{ii});
	% Convert to original scale:
	switch parameterNames{ii}
		case 'beta1'
			samples.(parameterNames{ii})		= z*sdy/sdx;
		case 'beta0'
			samples.(parameterNames{ii})		= z*sdy+muy-samples.beta1.*sdy*mux/sdx;
		case 'sigma'
			samples.(parameterNames{ii})	= z*sdy;
	end
end


%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

function pearson = genMCMCpearson(x,y)
%% Correlation Coefficient
chain_globals;

x		= [x y];

% Constants
[n,~]	= size(x);

%% Sampling
% MCMC Parameters
nburnin		= 500; % How Many Burn-in Samples?
nsamples	= 1000;  %How Many Recorded Samples?
nthin		= 1; % How Often is a Sample Recorded?

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n);

%% INTIALIZE THE CHAINS.
nChains		= nChainsDefault; % Number of chains to run.
initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).r		= 0;    % because data are standardized
	initsStruct(ii).mu		= zeros(1,2);        % because data are standardized
	initsStruct(ii).lambda	= ones(1,2);  % because data are standardized
end

modelname = which('pearson.txt');
% Use JAGS to Sample
samples = matjags( ...
	datastruct, ...
	modelname, ...
	initsStruct, ...
	'doparallel' , 1, ...
	'nchains', nChains,...
	'nburnin', nburnin,...
	'nsamples', nsamples, ...
	'thin', nthin, ...
	'monitorparams', {'r','mu','sigma'}, ...
	'savejagsoutput',0, ...
	'verbosity' ,0, ...
	'dic',0,...
	'cleanup',1, ...
	'workingdir' , 'tmpjags' );

r = [samples.r];
r = r(:);

mu		= [samples.mu];
[m,n,k] = size(mu);
mu		= reshape(mu,m*n,k);

sigma		= [samples.sigma];
[m,n,k]		= size(sigma);
sigma		= reshape(sigma,m*n,k);

pearson		= struct('r',r','mu',mu,'sigma',sigma);

function posteriorprediction(x,samples)
% POSTERIORPREDICTION(X,Y,B)
%
%

b0		= samples.beta0;
b1		= samples.beta1;
sigma	= samples.sigma;

% Posterior prediction:
% Specify x values for which predicted y's are needed:
xPostPred		= linspace(min(x),max(x),20);
% Define matrix for recording posterior predicted y values at each x value.
% One row per x value, with each row holding random predicted y values.
postSampSize	= length(b1);
yPostPred		= NaN(length(xPostPred),postSampSize);
% Define matrix for recording HDI limits of posterior predicted y values:
yHDIlim			= NaN(length(xPostPred),2);

% Generate posterior predicted y values.
% This gets only one y value, at each x, for each step in the chain.
for chainIdx = 1:postSampSize
	mu		= b0(chainIdx) + b1(chainIdx)*xPostPred;
	sd		= sigma(chainIdx);
	yPostPred(:,chainIdx) = mu+sd*randn(size(xPostPred));
end

for xIdx = 1:length(xPostPred)
	yHDIlim(xIdx,:) = hdimcmc(yPostPred(xIdx,:));
end

function writesimplemodel
% Placeholder function to generate JAGS model
str = ['model {\r\n',...
	'\tfor(i in 1:Ndata) {\r\n',...
	'\t\ty[i] ~ dnorm(mu[i],1/sigma)\r\n',...
	'\t\tmu[i] <- beta0 + beta1 * x[i]\r\n',...
	'\t}\r\n',...
	'\tbeta0 ~ dnorm( 0 , 1.0E-12 )\r\n',...
	'\tbeta1 ~ dnorm( 0 , 1.0E-12 )\r\n',...
	'\tsigma ~ dgamma(0.001,0.001 )\r\n',...
	'}\r\n',...
	];

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);

function writerobustmodel
% Placeholder function to generate JAGS model
str = ['model {\r\n',...
	'\tfor(i in 1:Ndata) {\r\n',...
	'\t\ty[i] ~ dt(mu[i],1/sigma^2,nu)\r\n',...
	'\t\tmu[i] <- beta0 + beta1 * x[i]\r\n',...
	'\t}\r\n',...
	'\tbeta0 ~ dnorm( 0 , 1.0E-12 )\r\n',...
	'\tbeta1 ~ dnorm( 0 , 1.0E-12 )\r\n',...
	'\tsigma ~ dunif(1.0E-3,1.0E+3 )\r\n',...
	'\tnu <- nuMinusOne+1\r\n',...
	'\tnuMinusOne ~ dexp(1/29.0)\r\n',...
	'}\r\n',...
	];

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);

function samples = extractchain(samples)
parameterNames	= fieldnames(samples); % get all parameter names
for ii			= 1:numel(parameterNames)
	n			= size(samples.(parameterNames{ii}));
	switch numel(n)
		case 2
			samples.(parameterNames{ii}) = samples.(parameterNames{ii})(:);
		case 3
			samples.(parameterNames{ii}) = reshape(samples.(parameterNames{ii}),n(1)*n(2),n(3));
			
	end
end