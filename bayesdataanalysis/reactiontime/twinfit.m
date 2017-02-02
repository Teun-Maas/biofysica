function samples = twinfit(t,RT,RT1,RT2,varargin)
% MCMC = TWINFIT(X,Y)
%
% Psychometric curve
%
% X = stimulus onset asynchrony
% Y = reaction time
%
% MCMC structure contains fields with MCMC samples:
% - lamba1: rate target 1
% - lambda2: rate distractor/target 2
% - mu: mean reaction time second stage
% - omega: window size
% - d: integration effect
% - k: warning effect
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
% Modified to Matlab code: Marc M. van Wanrooij
%
% See also MATJAGS, PLOTPOST
% TODO:
% - one gamma inferred for all subjects/conditions
% - Stan
% - determine guess rate from #alternatives
% - bias in guessing
% - JND (X>=0)



%% Initialization
% check for input orientation
% Y needs to be either:
% - a nx1 column vector with 0s and 1s
% - a nx2 column matrix with rate in colum 1 and #trials in column 2
if size(RT1,2)>2
	warning('Y should be a 1/2-column vector/matrix');
	RT1		= RT1';
	RT2		= RT2';
	RT		= RT';
	t		= t';
end

% model
task				= keyval('task',varargin,'FA'); % default: Focused Attention task

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,1000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,0000);
saveName		= keyval('saveName',varargin,'TWIN-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

% graphic flags
diagFlag		= keyval('showDiag',varargin,true); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
predFlag		= keyval('showPred',varargin,true); % show posterior predictive
centroidFlag	= keyval('showCentroid',varargin,'mode'); % mode, median, mean


%% Actual regression
samples				= genMCMC(t,RT,RT1,RT2,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic,'task',task);

%% Thin
if thinSteps>1
	samples = thinchain(samples,thinSteps);
end

%% MCMC diagnostics
if diagFlag
	parameterNames	= fieldnames(samples); % get all parameter names
	for parIdx			= 1:numel(parameterNames)
		n = size(samples.(parameterNames{parIdx}),3);
		for ii = 1:n
			figure
			a		= squeeze(samples.(parameterNames{parIdx})(:,:,ii));
			samp	= samples;
			samp.(parameterNames{parIdx}) = a;
			diagmcmc(samp,'parName',parameterNames{parIdx});
		end
	end
end

%% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

%% Posterior estimates
if postFlag
	figure;
	parameterNames	= fieldnames(samples); % get all parameter names
	n = numel(parameterNames);

	
	for parIdx			= 1:n
		a = samples.(parameterNames{parIdx});
		a = a(:);
		subplot(4,3,parIdx);
		cla
		plotpost(a);
		title(parameterNames{parIdx});
	end
end

%% Checks and display
if predFlag
	posteriorprediction(t,RT,RT1,RT2,samples,centroidFlag); % Posterior predictive
end

%% Sub-functions
function [samples,stats] = genMCMC(t,rt,rt1,rt2,varargin)
% B = GENMCMC(X,Y)
%
% Generate MCMC chains

%% initialization
% most is already defined in main function, so no defaults
numSavedSteps	= keyval('numSavedSteps',varargin);
thinSteps		= keyval('thinSteps',varargin);
burnInSteps		= keyval('burnInSteps',varargin);
saveName		= keyval('saveName',varargin);
nChains			= keyval('nChains',varargin);
runjagsMethod	= keyval('runjagsMethod',varargin);
dic				= keyval('dic',varargin);
modelname		= fullfile(pwd, 'model.txt');

%% Write the model
% first check parameters to be monitored

parameters		= {'omega','mu','delta','kappa','gamma','sigma','lambda1','lambda2'};

writemodel;

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.

RT = [rt; rt1; rt2];
T = [t'; repmat(t(end),size(rt1)); repmat(t(end),size(rt2))];
cond = [ones(size(rt)); repmat(2,size(rt1)); repmat(3,size(rt2))];
% Specify data, as a structure
dataStruct = struct('t',T,'rt',RT,'cond',cond,'nt',numel(T));

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lambda) are determined from dependent variables/
% y data
for ii = 1:nChains
	initsStruct(ii).mu			= 100; % ~mean
	initsStruct(ii).omega		= 300; % ~mean
	initsStruct(ii).kappa		= 50; % ~mean
	initsStruct(ii).delta		= 100; % ~mean
	initsStruct(ii).gamma		= 300; % ~mean
	initsStruct(ii).sigma		= 1; % ~mean
		initsStruct(ii).prompt1	= 100; % ~mean
	initsStruct(ii).prompt2		= 100; % ~mean

end

%% RUN THE CHAINS
% adaptSteps		= 500;			% Number of steps to 'tune' the samplers.
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethod,'parallel');
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

function posteriorprediction(t,rt,rt1,rt2,samples,centroidFlag)
% POSTERIORPREDICTION(X,Y,B)
%
%
% keyboard

%% Parameters

% 		samples.parameterNames
% 		gamma	= samples.gamma;
% delta	= samples.delta;
% kappa	= samples.kappa;
% omega	= samples.omega;
% mu		= samples.mu;
% sigma	= samples.sigma;


% Determine centroid
parameterNames	= fieldnames(samples); % get all parameter names
for parIdx			= 1:numel(parameterNames)
	str = [parameterNames{parIdx} ' = samples.(parameterNames{parIdx});'];
	eval(str);
	str = ['mu' parameterNames{parIdx} ' = ' centroidFlag '(samples.(parameterNames{parIdx}));'];
	eval(str);
end


%%
chainLength = length(omega);

% Graphics

% close all

figure;
subplot(131)
plot(t,rt,'ko','MarkerFaceColor','w');
hold on
% believable psychometruc curve
cVec		= floor(linspace(1,chainLength,30));
xWid		= max(t)-min(t);
xComb		= linspace(min(t)-0.1*xWid,max(t)+0.1*xWid,201);
xlim([min(t)-0.1*xWid,max(t)+0.1*xWid]);

for cIdx	= cVec
	ypred = twinfun(xComb,lambda1(cIdx),lambda2(cIdx),mu(cIdx),omega(cIdx),delta(cIdx),gamma(cIdx),kappa(cIdx));
	plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.7 .7 .7]);
end

function writemodel
% Placeholder function to generate JAGS model

%%
% str = ['# Specify the model for standardized data:\r\n',...
% 	'model {\r\n',...
% 	'\t# Bimodal\r\n',...
% 	'\tfor ( i in 1:nt ) {			\r\n',...
% 	'\t\t rt[i] ~ dnorm(rtmu[i],1/pow(sigma,2))	\r\n',...
% 	'\t\t rtmu[i] <- 1/lambda1+mu-delta*PI[i]-kappa*PW[i]		\r\n',...
% 	'\t\t# for every stimulus onset asynchrony t\r\n',...
% 	'\t\t# Probability of integration P(I)\r\n',...
% 	'\t\t PI[i] <- ifelse(x1[i],a[i],ifelse(x2[i],b[i],c[i]))		\r\n',...
% 	'\t\t x1[i] <- t[i]<(t[i]+omega) && (t[i]+omega)<0	\r\n',...
% 	'\t\t x2[i] <- t[i]<=0 && 0<=(t[i]+omega)		\r\n',...
% 	'\t\t# x3[i] <- 0<t[i] && t[i]<(t[i]+omega)		\r\n',...
% 	'\t\t a[i] <- lambda1/(lambda1+lambda2)*exp(lambda2*t[i])*(-1+exp(lambda2*omega))				\r\n',...
% 	'\t\t b[i] <- 1/(lambda1+lambda2)*(lambda2*(1-exp(-lambda1*(omega+t[i])))+lambda1*(1-exp(lambda2*t[i]))) \r\n',...
% 	'\t\t c[i] <- lambda2/(lambda1+lambda2)*(exp(-lambda1*t[i])-exp(-lambda1*(omega+t[i])))			\r\n',...
% 	'\t\t# Probability of warning P(W)\r\n',...
% 	'\t\t PW[i] <- ifelse(x4[i],d[i],e[i])		\r\n',...
% 	'\t\t x4[i] <- (t[i]+gamma)<0	\r\n',...
% 	'\t\t d[i] <- 1-(lambda1/(lambda1+lambda2)*exp(lambda2*(t[i]+gamma)))			\r\n',...
% 	'\t\t e[i] <- (lambda2/(lambda1+lambda2)*exp(-lambda1*(t[i]+gamma)))	\r\n',...
% 	'\t}\r\n',...
% 	'\t# Bimodal priors\r\n',...
% 	'\tmu ~ dnorm(0 , 1.0E-12)\r\n',...
% 	'\tsigma ~ dgamma(1,1)\r\n',...
% 	'\tomega ~ dgamma(1,1)\r\n',...
% 	'\tdelta ~ dgamma(1,1)\r\n',...
% 	'\tkappa ~ dgamma(1,1)\r\n',...
% 	'\tgamma ~ dgamma(1,1)\r\n',...
% 	'\tlambda1 ~ dgamma(1,1)\r\n',...
% 	'\tlambda2 ~ dgamma(1,1)\r\n',...
% 	'\t}\r\n',...
% 	];


str = ['# Specify the model for standardized data:\r\n',...
	'model {\r\n',...
	'\t# Bimodal\r\n',...
	'\tfor ( i in 1:nt ) {			\r\n',...
	'\t\t rt[i] ~ dnorm(rtmu[i],1/pow(sigma,2))	\r\n',...
	'\t\t rtmu[i] <- ifelse(cond[i]==1,rtmubi[i],ifelse(cond[i]==2,rtmu1[i],rtmu2[i]))		\r\n',...
	'\t\t rtmubi[i] <- 1/lambda1+mu-delta*PI[i]-kappa*PW[i]		\r\n',...
	'\t\t# rtmubi[i] <- 1/lambda1+mu-delta*PI[i]-kappa*PW[i]		\r\n',...
	'\t\t rtmu1[i] <- 1/lambda1+mu		\r\n',...
	'\t\t rtmu2[i] <- 1/lambda2+mu		\r\n',...
	'\t\t# for every stimulus onset asynchrony t\r\n',...
	'\t\t# Probability of integration P(I)\r\n',...
	'\t\t PI[i] <- ifelse(x1[i],a[i],ifelse(x2[i],b[i],c[i]))		\r\n',...
	'\t\t x1[i] <- t[i]<(t[i]+omega) && (t[i]+omega)<0	\r\n',...
	'\t\t x2[i] <- t[i]<=0 && 0<=(t[i]+omega)		\r\n',...
	'\t\t# x3[i] <- 0<t[i] && t[i]<(t[i]+omega)		\r\n',...
	'\t\t a[i] <- lambda1/(lambda1+lambda2)*exp(lambda2*t[i])*(-1+exp(lambda2*omega))				\r\n',...
	'\t\t b[i] <- 1/(lambda1+lambda2)*(lambda2*(1-exp(-lambda1*(omega+t[i])))+lambda1*(1-exp(lambda2*t[i]))) \r\n',...
	'\t\t c[i] <- lambda2/(lambda1+lambda2)*(exp(-lambda1*t[i])-exp(-lambda1*(omega+t[i])))			\r\n',...
	'\t\t# Probability of warning P(W)\r\n',...
	'\t\t PW[i] <- ifelse(x4[i],d[i],e[i])		\r\n',...
	'\t\t x4[i] <- (t[i]+gamma)<0	\r\n',...
	'\t\t d[i] <- 1-(lambda1/(lambda1+lambda2)*exp(lambda2*(t[i]+gamma)))			\r\n',...
	'\t\t e[i] <- (lambda2/(lambda1+lambda2)*exp(-lambda1*(t[i]+gamma)))	\r\n',...
	'\t}\r\n',...
	'\t# Bimodal priors\r\n',...
	'\tmu ~ dgamma(1,1)\r\n',...
	'\tsigma ~ dunif(1.0E-3,1)\r\n',...
	'\tomega ~ dunif(0,600)\r\n',...
	'\tdelta ~ dunif(0,600)\r\n',...
	'\tkappa ~ dgamma(1,1)\r\n',...
	'\tgamma ~ dgamma(1,1)\r\n',...
	'\tlambda1 <- 1/prompt1\r\n',...
	'\tlambda2 <- 1/prompt2\r\n',...
	'\tprompt1 ~ dunif(5,150)\r\n',...
	'\tprompt2 ~ dunif(5,150)\r\n',...	
	'\t}\r\n',...
	];
% dunif(1.0E-3,1.0E+3 )

% 	parstart = [lv		la		m	o	d	ga	k;
% 				1/150	1/150	0	0	0	0	0;
% 				1/5		1/5		Inf 600 Inf Inf Inf];
			
% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);
