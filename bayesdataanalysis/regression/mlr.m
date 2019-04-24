function [esamples,stats] = mlr(y,x,varargin)
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
	x1 = -90:10:90;
	x2 = 40:1:70;
	[x1,x2] = meshgrid(x1,x2);
	x1 = x1(:);
	x2 = x2(:);
	x = [x1 x2];
	
% 	x = x';
	y = 0.9*x1+0.4*x2+10*randn(size(x1));
	
		whos x1 x y

end

%% Initialization
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,500); 
thinSteps		= keyval('thinSteps',varargin,1);
burnInSteps		= keyval('burnInSteps',varargin,100);
saveName		= keyval('saveName',varargin,'HierMultLinRegressData-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
type			= keyval('type',varargin,'simple'); % simple or robust
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
diagFlag		= keyval('diag',varargin,false);
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
% keyboard
%% Extract chain values:
esamples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

%
% beta1 = esamples.beta(:,1);
% beta2 = esamples.beta(:,2);
% 
% close all
% subplot(121)
% plotpost(beta1,'xlim',[0 1.2]);
% 
% subplot(122)
% plotpost(beta2,'xlim',[0 1.2]);

%% Pearson correlation
% pearson		= genMCMCpearson(x,y);
% 
% % samples		= addfields(samples,pearson);
% samples.pearson = pearson;

%% Checks and display
% posteriorprediction(x,samples);

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
dic				= keyval('dic',varargin);

% modelname = which('regress.txt');
% if ~exist(modelname,'file')
% if strcmpi(type,'simple')
% 	writesimplemodel;
% 	modelname = fullfile(pwd, 'model.txt');
% 	parameters		= {'beta0','beta1','sigma'};		% The parameter(s) to be monitored.
% elseif strcmpi(type,'robust')
% 	writerobustmodel;
% 	modelname = fullfile(pwd, 'model.txt');
% 	parameters		= {'beta0','beta1','sigma','nu'};		% The parameter(s) to be monitored.
% 	
% end

% modelname = fullfile(pwd, 'mlr.txt');
modelname = which('mlr.txt')
parameters = {'beta0' ,  'beta' ,  'sigma', 'zbeta0' , 'zbeta' , 'zsigma', 'nu'};
% end

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
% [zx,mux,sdx]		= zscore(x);
% [zy,muy,sdy]		= zscore(y);
% nData				= size(x,1);

[Ntotal,Nx] = size(x);

% Specify data, as a structure
dataStruct = struct('x',x,'y',y,'Nx',Nx','Ntotal',Ntotal);

%   dataList = list(
%     x = x ,
%     y = y ,
%     Nx = dim(x)[2] ,
%     Ntotal = dim(x)[1]
%   )
%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
parameters = {'beta0' ,  'beta' ,  'sigma', 'zbeta0' , 'zbeta' , 'zsigma', 'nu'};

for ii = 1:nChains
	initsStruct(ii).zbeta0		= 0;    % because data are standardized
% 	initsStruct(ii).beta1		= 0;        % because data are standardized
% 	initsStruct(ii).sigma		= (1-r(2)^2);  % because data are standardized
% 		initsStruct(ii).nu		= 100;        % because data are standardized
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



%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

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



