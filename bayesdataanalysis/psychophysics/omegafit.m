function samples = omegafit(x,y,s,varargin)
% MCMC = OMEGAFIT(X,Y,S)
%
% Psychometric curve with ordinal data
%
% X = signal level
% Y = 0 or 1 (incorrect/correct) or Y = [rate ntrials]
% S = subject ID (or condition)
%
% MCMC structure contains fields with MCMC samples:
% - beta0: offset
% - beta1: slope
% - guess: guess rate
%
% MCMC = REGJAGS(Y,X,'NAME',VALUE)
% Additional name-value pair inputs include:
% - 'type' - 'logistic' (default)
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
%
% See also
% Kuss, M., J?kel, F., & Wichmann, F. A. (2005). Bayesian inference for
% psychometric functions. Journal of Vision, 5(5), 478?92.
% doi:10.1167/5.5.8
%
% Modified to Matlab code: Marc M. van Wanrooij
%
% See also MATJAGS, PLOTPOST
% TODO:
% - one gamma inferred for all subjects/conditions
% - Stan
% - determine guess rate from #alternatives
% - bias in guessing
% - JND (X>=0)
% predFlag		= keyval('showPred',varargin,false); % show posterior predictive



%% Initialization
% x = linspace(-10,20,100);
% p = normpdf(x,(1+6)/2,(1+6)/2);
% close all
% plot(x,p);
% xlim([1 6]);
% ylim([0 max(p)*1.1]);
% return
%% if no data, use DBDA data
if nargin<1

	
	
% 	  mu[i] <- alpha*exp(-(zx[i]-zgmu)^2/(2*zgsigma^2))
fun		= @(beta,x)( beta(1).*exp( -1/2 .* (x-beta(2)).^2 ./ beta(3)^2 ) );
x		= linspace(-25,10,20);
nrep	= 2;
x		= repmat(x,1,nrep);

alpha	= 9.3;
theta	= -11;
gsigma	= 6;
sigma	= 0.3;

beta	= [alpha theta gsigma];
y		= fun(beta,x)+sigma*randn(size(x));



scores		= 1:10;
nLevel		= max(scores);
threshold	= 1.5:(nLevel-0.5);

score		= y;
nthresh = numel(threshold);
sel = y<=threshold(1);
N(1) = sum(sel);
score(sel) = 1;
for ii = 2:nthresh
	sel = y<=threshold(ii) & y>threshold(ii-1);
	N(ii) = sum(sel);
	score(sel) = ii;

end
sel = y>threshold(nthresh);
N(nLevel) = sum(sel);
score(sel) = nLevel;

y = score;
end
% return
% check for input orientation
% Y needs to be either:
% - a nx1 column vector with 0s and 1s
% - a nx2 column matrix with rate in colum 1 and #trials in column 2
% if size(y,2)>2
% 	warning('Y should be a 1/2-column vector/matrix');
% 	y = y';
% 	x = x';
% 	s = s';
% end
% if size(y,1)~=size(x,1)
% 	warning('Y and X should have same orientations');
% 	x = x';
% end
% if size(y,1)~=size(x,1)
% 	error('X and Y should have same number of rows');
% end
% if isempty(s)
% 	s = ones(size(y,1),1);
% end
%% S should contain consecutive values, no missing values
[~,~,s] = unique(s);

% % model
% if size(y,2)==2
% 	respDist = 'binomial';
% elseif size(y,2)==1
% 	respDist = 'bernouilli';
% end


% fun				= keyval('function',varargin,@logisticfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
fig				= keyval('figure',varargin,200); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun

% gamma			= keyval('gamma',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 gamma for all subjects/conditions
% lambda			= keyval('lambda',varargin,'infer');% 0.001 = low lapse rate, 'infer' = estimate from data, 'single' = infer 1 lambda for all subjects/conditions
nLevels			= keyval('nLevels',varargin,max(y)); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,5000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,1000);
saveName		= keyval('saveName',varargin,'Hier-OrdinalCurve-Jags-');
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

samples = jags_omegafit(x,y,s,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic,'nLevels',nLevels);


%% Thin
if thinSteps>1
	samples = thinchain(samples,thinSteps);
end

%% MCMC diagnostics
if diagFlag
	diagmcmc(samples)
end

% keyboard
%% Extract chain values:
samples = extractchain(samples); % from multiple -dimensiorn matrix to 1- or 2-D




function [samples,stats] = jags_omegafit(x,y,s,varargin)
% SAMPLES = JAGS_PSIFIT(X,Y)
%
% See also PSIFIT

%% Psychometric function
modelname = '/Users/marcw/Documents/GitLab/biofysica/bayesdataanalysis/psychophysics/jags_omega_model.txt';
%% MCMC parameters
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-OmegaCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
nLevels			= keyval('nLevels',varargin,max(y)); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
%% Numbers (amount of data and groups)
[Nx,Ntotal]		= size(x)

Nsubj			= numel(unique(s));
%% Data
% parameters			= {'thresh','mu','sigma'};
parameters			= {'alpha','sigma','thresh','theta','omega','mualpha','mutheta','muomega','musigma'};
% parameters			= {'alpha','sigma','theta','omega','mualpha','mutheta','muomega'};

thresh				= NaN(nLevels-1,1);
thresh(1)			= 1 + 0.5;
thresh(nLevels-1)	= nLevels-1 + 0.5;
% thresh = 1.5:(nLevels-1 + 0.5);
dataStruct			= struct('x',x,'y',y,'s',s','Nsubj',Nsubj,'Ntotal',Ntotal,'nYlevels',nLevels,'thresh',thresh); % data structure
initsStruct			= initialize(nChains,x,y,s); % initial parameter values

%% parallel?
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end

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



function initsStruct	= initialize(nChains,x,y,s)
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lapse rate) are determined from dependent variables/
% y data

[mx,idx]	= max(y);
S			= assemble(y,x,'fun',@std);
Nsubj		= max(unique(s));
for ii = 1:nChains
	
% 	initsStruct(ii).omega			= repmat(std(x),Nsubj,1);  % standardized std
% 	initsStruct(ii).theta			= repmat(x(idx),Nsubj,1);  % standardized std
	initsStruct(ii).mualpha			= mx;  % standardized std
if mean(S)>0.0001
	initsStruct(ii).musigma			= mean(S)/Nsubj;  % standardized std
else
	initsStruct(ii).musigma			= 0.0001;
end
	initsStruct(ii).mutheta			= x(idx);  % standardized std

end



