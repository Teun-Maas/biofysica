function [samples,stats] = jags_psifit(x,y,s,varargin)
% SAMPLES = JAGS_PSIFIT(X,Y)
%
% See also PSIFIT

%% Psychometric function
fun				= keyval('function',varargin,@logisticfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
gamma			= keyval('gamma',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 gamma for all subjects/conditions
lambda			= keyval('lambda',varargin,'infer');% 0.001 = low lapse rate, 'infer' = estimate from data, 'single' = infer 1 lambda for all subjects/conditions
modelname		= fullfile(pwd, 'jags_psi_model.txt');

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
writemodel(gamma,lambda,Ngroups,fun,respDist,modelname);

%% Data
parameters		= checkwhich(gamma,lambda,Ngroups);
[zx,mux,sdx]	= standardize(x,fun);
dataStruct		= struct('x',zx,'y',y,'s',s,'Ntotal',Ntotal,'Nsubj',Ngroups); % data structure
initsStruct		= initialize(Ngroups,nChains,gamma,lambda,fun); % initial parameter values

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


samples				= unzscore(samples,mux,sdx);

%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

function [zx,mux,sdx]	= standardize(x,fun)

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
[zx,mux,sdx]		= zscore(x);
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	% Weibull functions are defined for x>0
	mux = 0;
	sdx = 1;
	zx	= (x-mux)/sdx;
end

function parameters		= checkwhich(gamma,lambda,ngroups)

if isnumeric(gamma) && isnumeric(lambda)&& ngroups>1
	parameters		= {'theta','omega','mutheta','muomega'};
elseif  isnumeric(gamma) && isnumeric(lambda)&& ngroups==1
		parameters		= {'theta','omega'};
elseif isnumeric(gamma) && ~isnumeric(lambda) && ngroups>1
	parameters		= {'theta','omega','lambda','mulambda','mutheta','muomega'};
elseif isnumeric(gamma) && ~isnumeric(lambda) && ngroups==1
	parameters		= {'theta','omega','lambda'};
elseif ~isnumeric(gamma) && isnumeric(lambda) && ngroups>1
	parameters		= {'theta','omega','gamma','mugamma','mutheta','muomega'};
elseif ~isnumeric(gamma) && isnumeric(lambda) && ngroups==1
	parameters		= {'theta','omega','gamma'};
elseif ~isnumeric(gamma) && ~isnumeric(lambda) && ngroups>1
	parameters		= {'theta','omega','gamma','lambda','mutheta','muomega','mugamma','mulambda'};
else
	parameters		= {'theta','omega','gamma','lambda'};
end

function initsStruct	= initialize(Nsubj,nChains,gamma,lambda,fun)
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lapse rate) are determined from dependent variables/
% y data
for ii = 1:nChains
	if strcmp(func2str(fun),'weibullfun')
	initsStruct(ii).theta			= ones(Nsubj,1); % standardized mean
	initsStruct(ii).omega			= ones(Nsubj,1);  % standardized std
	else
	initsStruct(ii).theta			= zeros(Nsubj,1); % standardized mean
	initsStruct(ii).omega			= ones(Nsubj,1);  % standardized std
	end
	if ~isnumeric(gamma)
		initsStruct(ii).gamma		= repmat(0.01,Nsubj,1); % no guessing
	end
	if ~isnumeric(lambda)
		initsStruct(ii).lambda		= repmat(0.01,Nsubj,1); % no lapses
	end
end

function samples		= unzscore(samples,mux,sdx)
parameterNames		= fieldnames(samples); % get all parameter names
for ii				= 1:numel(parameterNames)
	z			= samples.(parameterNames{ii});
	% Convert to original scale:
	switch parameterNames{ii}
		case 'theta'
			samples.(parameterNames{ii})		= z*sdx+mux;
		case 'mutheta'
			samples.(parameterNames{ii})		= z*sdx+mux;
		case 'omega'
			samples.(parameterNames{ii})		= z*sdx;
		case 'muomega'
			samples.(parameterNames{ii})		= z*sdx;
	end
end

function writemodel(gamma,lambda,ngroups,fun,respDist,modelname)
% Placeholder function to generate JAGS model


%% Psychometric function
psifun = func2str(fun);

switch psifun
	case 'logisticfun'
		fun = 'ilogit'; % 1/(1+exp(x)), 
		par = '(2*log(1/alpha-1))/omega[s[i]] * (x[i]-theta[s[i]])';% eqn. 7 Kuss et al 2005
	case 'probitfun'
		fun = 'phi';
		par = '( (probit(1-alpha)-probit(alpha))/omega[s[i]] ) * (x[i]-theta[s[i]])';
	case 'gumbelfun'
		fun = 'icloglog';
		par = '( (log(-log(alpha))-log(-log(1-alpha)))/omega[s[i]] * (x[i]-theta[s[i]]) + log(-log(0.5)) )';
	case 'revgumbelfun'
		fun = '1-icloglog';
		par = '( (log(-log(1-alpha))-log(-log(alpha)))/omega[s[i]] * (x[i]-theta[s[i]]) + log(-log(0.5)) )';
	case 'weibullfun'
		fun = 'icloglog'; % 1-exp(-exp(x))
		par = '( 2*omega[s[i]]*theta[s[i]]/log(2)  * ( log(x[i]) - log(theta[s[i]]) ) + log(log(2)) )'; % eqn. 11 Kuss et al 2005
	case 'revweibullfun'
		fun = '1-icloglog';
		par = '( -2*omega[s[i]]*theta[s[i]]/log(2)  * ( log(x[i]) - log(theta[s[i]]) ) + log(log(2)) )';
end

%% Guess rate gamma
if isnumeric(gamma)
	gammastrmu	= 'gamma';
	gammastr1	= 	[];
	gammastr2	= 	['\t gamma\t\t\t<- ' num2str(gamma) '  \r\n'];
elseif strcmp(gamma,'infer')
	gammastrmu = 'gamma[s[i]]';
	gammastr1	=	'\t\t gamma[j]\t~ dbeta(agamma,bgamma)  \r\n';
	gammastr2	= [	'\t# guess rate / gamma hyperprior \r\n',...
		'\t agamma\t\t\t<- mugamma*kappagamma \r\n',...
		'\t bgamma\t\t\t<- (1.0-mugamma) * kappagamma  \r\n',...
		'\t mugamma\t\t~ dbeta(2.0,2.0)  \r\n',...
		'\t kappagamma\t\t~ dgamma(Skappa,Rkappa)  \r\n'];
	if ngroups==1
			gammastr2	= [	'\t# guess rate / gamma hyperprior \r\n',...
		'\t agamma\t\t\t<- 1.0 \r\n',...
		'\t bgamma\t\t\t<- 1.0  \r\n'];
	end
end

%% Lapse rate lambda
if isnumeric(lambda)
	lambdastrmu	= 'lambda';
	lambdastr1	= 	[];
	lambdastr2	= 	['\t lambda\t\t\t<- ' num2str(lambda) '  \r\n'];
	
elseif strcmp(lambda,'infer')
	lambdastrmu = 'lambda[s[i]]';
	lambdastr1	=	'\t\t lambda[j]\t~ dbeta(alambda,blambda)  \r\n';
	lambdastr2	= [	'\t# lapse rate/lambda hyperprior \r\n',...
		'\t alambda\t\t<- mulambda*kappalambda \r\n',...
		'\t blambda\t\t<- (1.0-mulambda) * kappalambda  \r\n',...
		'\t mulambda\t\t~ dbeta(2.0,2.0)  \r\n',...
		'\t kappalambda\t~ dgamma(Skappa,Rkappa)  \r\n'];
	if ngroups==1
		lambdastr2	= [	'\t# lapse rate/lambda hyperprior \r\n',...
			'\t alambda\t\t<- 1.0 \r\n',...
			'\t blambda\t\t<- 1.0  \r\n'];
	end
end

if isempty(lambdastr1) && isempty(gammastr1)
	hypergammastr = [];
else
	hypergammastr = [	'\t# constants for all dgamma hyperprior \r\n',...
		'\t Skappa\t\t\t<- pow(10,2)/pow(10,2) # constant \r\n',...
		'\t Rkappa\t\t\t<- 10/pow(10,2) # constant \r\n',...
		];
end
%% Mean rate
mustr		=	['\t\t mu[i]\t<-  (1-' lambdastrmu ') * ((1-' gammastrmu ') * (' fun '( ' par ' )) + ' gammastrmu ' ) + (' gammastrmu '*' lambdastrmu ')\r\n'];

%% Response distribution
switch respDist
	case 'bernouilli'
		distr		= '\t\t y[i] ~ dbern(mu[i])\r\n';
	case 'binomial'
		distr		= '\t\t y[i,1] ~ dbin(mu[i],y[i,2])\r\n';
end

%% Hyperpriors
if ngroups>1
	if strcmp(psifun,'weibullfun')
hyperprior = ['\t# threshold/theta hyperprior \r\n',...
			'\t mutheta\t\t\t<- 1.0E-3\r\n',...
			'\t sigmatheta\t\t<- 1000\r\n',...
			'\t# width/omega hyperprior \r\n',...
			'\t somega\t\t\t<- 1.0E-3 # shape \r\n',...
			'\t romega\t\t\t<- 1000 # rate \r\n',...
		'\t muomega\t\t\t\t<-somega\r\n',...
		'\t sigmaomega\t\t\t\t<-romega \r\n'];
	else
		hyperprior = ['\t# threshold/theta hyperprior \r\n',...
		'\t mutheta\t\t\t~ dnorm(0,1/(10)^2)\r\n',...
		'\t sigmatheta\t\t~ dunif(1.0E-3,1000)\r\n',...
		'\t# width/omega hyperprior \r\n',...
		'\t somega\t\t\t<- pow(muomega,2)/pow(sigmaomega,2) # shape \r\n',...
		'\t romega\t\t\t<- muomega/pow(sigmaomega,2) # rate \r\n',...
		'\t muomega\t\t\t\t~ dgamma(1,0.25)\r\n',...
		'\t sigmaomega\t\t\t\t~ dgamma(1,0.5) \r\n'];
	end
else
	if strcmp(psifun,'weibullfun')
		hyperprior = ['\t# threshold/theta hyperprior \r\n',...
			'\t mutheta\t\t\t<- 1.0E-3\r\n',...
			'\t sigmatheta\t\t<- 1000\r\n',...
			'\t# width/omega hyperprior \r\n',...
			'\t somega\t\t\t<- 1.0E-3 # shape \r\n',...
			'\t romega\t\t\t<- 1000 # rate \r\n'];
	else
		hyperprior = ['\t# threshold/theta hyperprior \r\n',...
			'\t mutheta\t\t\t<- 0\r\n',...
			'\t sigmatheta\t\t<- 1/(10)^2\r\n',...
			'\t# width/omega hyperprior \r\n',...
			'\t somega\t\t\t<- 1.0 # shape \r\n',...
			'\t romega\t\t\t<- 0.25 # rate \r\n'];
	end
end

%% priors
omegapriorstr =	'\t\t omega[j]\t~ dgamma(somega,romega)  \r\n';
thetapriorstr =	'\t\t theta[j]\t~ dnorm(mutheta,sigmatheta)\r\n';
if strcmp(psifun,'weibullfun')
	thetapriorstr =	'\t\t theta[j]\t~ dunif(mutheta,sigmatheta)\r\n';
	omegapriorstr =	'\t\t omega[j]\t~ dunif(somega,romega)  \r\n';
end


%% Model
str = ['# Specify the model for standardized data:\r\n',...
	'model {\r\n',...
	'\tfor ( i in 1:Ntotal ) {\r\n',...
	'\t\t# Bernouilli distributed trials\r\n',...
	distr,...
	mustr,...
	'\t}\r\n',...
	'\tfor ( j in 1:Nsubj ) {\r\n',...
	'\t\t# Priors vague on standardized scale:\r\n',...
	gammastr1,...
	lambdastr1,...
	omegapriorstr,...
	thetapriorstr,...
	'\t}\r\n',...
	'\t# constant/default width at F-1(ALPHA) and F-1(1-ALPHA) \r\n',...
	'\t alpha\t\t\t<- 0.1 \r\n',...
	gammastr2,...
	lambdastr2,...
	hypergammastr,...
	hyperprior,...
	'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen(modelname,'w');
fprintf(fid,str);
fclose(fid);

