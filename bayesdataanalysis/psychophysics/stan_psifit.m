function [samples,stats] = stan_psifit(x,y,s,varargin)
% SAMPLES = JAGS_PSIFIT(X,Y)
%
% See also PSIFIT

%% Psychometric function
fun				= keyval('function',varargin,@logisticfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
gamma			= keyval('gamma',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 gamma for all subjects/conditions
lambda			= keyval('lambda',varargin,'infer');% 0.001 = low lapse rate, 'infer' = estimate from data, 'single' = infer 1 lambda for all subjects/conditions
modelname		= fullfile(pwd, 'stan_psi_model.stan');

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
model = writemodel(gamma,lambda,Ngroups,fun,respDist,modelname);

%% Data
parameters		= checkwhich(gamma,lambda,Ngroups);
[zx,mux,sdx]	= standardize(x,fun);
% dataStruct		= struct('x',zx,'y',y,'s',s,'N',Ntotal,'Nsubj',Ngroups); % data structure
dataStruct		= struct('x',zx,'y',y,'N',Ntotal,'alpha',0.1); % data structure
if isnumeric(gamma)
 dataStruct.('g')=gamma ; % data structure
end
if isnumeric(lambda)
 dataStruct.('l')=lambda ; % data structure
end
initsStruct		= initialize(Ngroups,nChains,gamma,lambda); % initial parameter values



%% parallel?
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end

%% MCMC
fprintf( 'Running STAN...\n' );
% [samples, stats] = matjags( ...
% 	dataStruct, ...                     % Observed data
% 	modelname, ...    % File that contains model definition
% 	initsStruct, ...                          % Initial values for latent variables
% 	'doparallel' , doparallel, ...      % Parallelization flag
% 	'nchains', nChains,...              % Number of MCMC chains
% 	'nburnin', burnInSteps,...              % Number of burnin steps
% 	'nsamples', nIter, ...           % Number of samples to extract
% 	'thin', thinSteps, ...                      % Thinning parameter
% 	'dic',dic, ...                       % Do the DIC?
% 	'monitorparams', parameters, ...     % List of latent variables to monitor
% 	'savejagsoutput',0, ...          % Save command line output produced by JAGS?
% 	'verbosity',0, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
% 	'cleanup',1);                    % clean up of temporary files?

% fit = stan( ...
% 	'model_code',model,...		% File that contains model definition
% 	'data',dataStruct,...		% Observed data
% 	'iter',nIter,...			% Number of samples to extract
% 	'chains',nChains,...		% Number of MCMC chains
% 	'warmup',burnInSteps,...	% Number of burnin steps
% 	'thin', thinSteps, ...      % Thinning parameter
% 	'model_name','stan_psi_model',...
% 	'file_overwrite',true);


fit = stan( ...
	'model_code',model,...		% File that contains model definition
	'data',dataStruct,...		% Observed data
	'init',initsStruct,...		% Observed data
	'iter',1000,...			% Number of samples to extract
	'chains',4,...		% Number of MCMC chains
	'warmup',500,...	% Number of burnin steps
	'thin', 1, ...      % Thinning parameter
	'model_name','stan_psi_model',...
	'file_overwrite',true,...
	'working_dir',pwd);

% fit = stan('model_code',model,'data',dataStruct,'model_name','stan_psi_model','file_overwrite',true);
% fit = stan('fit',fit,'data',dataStruct,'iter',10000,'chains',4);
fit.block();

samples				= stan2jagsmat(fit);

samples				= unzscore(samples,mux,sdx);

%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

stats = [];

function [zx,mux,sdx]	= standardize(x,fun)

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
[zx,mux,sdx]		= zscore(x);
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	% Weibull functions are defined for x>0
	mux = min(x);
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

function initsStruct	= initialize(Nsubj,nChains,gamma,lambda)
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lapse rate) are determined from dependent variables/
% y data
for ii = 1:nChains
	initsStruct(ii).theta			= zeros(Nsubj,1); % standardized mean
	initsStruct(ii).omega			= ones(Nsubj,1);  % standardized std
	if ~isnumeric(gamma)
		initsStruct(ii).g		= zeros(Nsubj,1); % no guessing
	end
	if ~isnumeric(lambda)
		initsStruct(ii).l		= zeros(Nsubj,1); % no lapses
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

function model = writemodel(gamma,lambda,ngroups,fun,respDist,modelname)
%% Psychometric function


%%


%% Psychometric function
fun = @logisticfun
gamma = 'infer';
lambda = 0;
respDist = 'bernouilli';
ngroups = 1;

fun = func2str(fun);
switch fun
	case "logisticfun"
		fun = "inv_logit";
		par = "(2*log(1/alpha-1))/omega * (x-theta)";
	case "probitfun"
		fun = "phi";
		par = "( (probit(1-alpha)-probit(alpha))/omega[s[i]] ) * (x[i]-theta[s[i]])";
	case "gumbelfun"
		fun = "icloglog";
		par = "( (log(-log(alpha))-log(-log(1-alpha)))/omega[s[i]] * (x[i]-theta[s[i]]) + log(-log(0.5)) )";
	case "revgumbelfun"
		fun = "1-icloglog";
		par = "( (log(-log(1-alpha))-log(-log(alpha)))/omega[s[i]] * (x[i]-theta[s[i]]) + log(-log(0.5)) )";
	case "weibullfun"
		fun = "icloglog";
		par = "( (2*omega[s[i]]*theta[s[i]])/log(2) * (log(x[i])-log(theta[s[i]])) + log(log(2)) )";
	case "revweibullfun"
		fun = "1-icloglog";
		par = "( (-2*omega[s[i]]*theta[s[i]])/log(2) * (log(x[i])-log(theta[s[i]])) + log(log(2)) )";
end

% Guess rate gamma
if isnumeric(gamma)
	gammastrmu	= "g";
	gammastr1	= 	[];
	gammastr2	= 	" gamma   <- "+num2str(gamma);
	fixgammastr	= 	"  real<lower=0,upper=1> g;";
	infgammastr	= 	[];
elseif strcmp(gamma,"infer")
	fixgammastr	= 	[];
	gammastrmu	= "g";
	gammastr1	=	" g ~ beta(1.0,1.0);";
	infgammastr	= 	"  real<lower=0,upper=1> g;";

	gammastr2	= [	" # guess rate / gamma hyperprior  "
		" agamma   <- mugamma*kappagamma  "
		" bgamma   <- (1.0-mugamma) * kappagamma   "
		" mugamma  ~ dbeta(2.0,2.0)   "
		" kappagamma  ~ dgamma(Skappa,Rkappa)   "];
	if ngroups==1
			gammastr2	= [	" # guess rate / gamma hyperprior  "
		" agamma   <- 1.0  "
		" bgamma   <- 1.0   "];
	end
end

% Lapse rate lambda
if isnumeric(lambda)
	lambdastrmu	= "l";
	lambdastr1	= 	[];
	lambdastr2	= 	" l   <- " + num2str(lambda);
	fixlambdastr	= 	"  real<lower=0,upper=1> l;";
elseif strcmp(lambda,"infer")
	fixlambdastr	= 	[];
	lambdastrmu = "lambda[s[i]]";
	lambdastr1	=	"  lambda[j] ~ dbeta(alambda,blambda)   ";
	lambdastr2	= [	" # lapse rate/lambda hyperprior  "
		" alambda  <- mulambda*kappalambda  "
		" blambda  <- (1.0-mulambda) * kappalambda   "
		" mulambda  ~ dbeta(2.0,2.0)   "
		" kappalambda ~ dgamma(Skappa,Rkappa)   "];
	if ngroups==1
		lambdastr2	= [	" # lapse rate/lambda hyperprior  "
			" alambda  <- 1.0  "
			" blambda  <- 1.0   "];
	end
end

if isempty(lambdastr1) && isempty(gammastr1)
	hypergammastr = [];
else
	hypergammastr = [	" # constants for all dgamma hyperprior  "
		" Skappa   <- pow(10,2)/pow(10,2) # constant  "
		" Rkappa   <- 10/pow(10,2) # constant  "
		];
end
% Mean rate
mustr		=	" rate = (1-"+lambdastrmu+") * ((1-"+gammastrmu+") * ("+fun+"( "+par+" )) + "+gammastrmu+" ) + ("+gammastrmu+"*"+lambdastrmu+"); ";
% mustr		=	" rate = "+fun+"("+par+");";

% Response distribution
switch respDist
	case 'bernouilli'
		distr		= " y ~ bernoulli(rate);";
	case 'binomial'
		distr		= " y ~ binomial(rate,k);";
end

% Hyperpriors
if ngroups>1
hyperprior = [" # threshold/theta hyperprior  "
	" mutheta   ~ dnorm(0,1/(10)^2) "
	" sigmatheta  ~ dunif(1.0E-3,1000) "
	" # width/omega hyperprior  "
	" somega   <- pow(muomega,2)/pow(sigmaomega,2) # shape  "
	" romega   <- muomega/pow(sigmaomega,2) # rate  "
	" muomega    ~ dgamma(1,0.25) "
	" sigmaomega    ~ dgamma(1,0.5)  "];
else
hyperprior = [" # threshold/theta hyperprior  "
	" mutheta   <- 0 "
	" sigmatheta  <- 1/(10)^2 "
	" # width/omega hyperprior  "
	" somega   <- 1.0 # shape  "
	" romega   <- 0.25 # rate  "];
end
% Model
str = [
	"# Specify the model for standardized data: "
	"model { "
	" for ( i in 1:Ntotal ) {"
	"  # Bernouilli distributed trials "
	distr
	mustr
	" } "
	" for ( j in 1:Nsubj ) { "
	"  # Priors vague on standardized scale: "
	gammastr1
	lambdastr1
	"  omega[j] ~ dgamma(somega,romega)   "
	"  theta[j] ~ dnorm(mutheta,sigmatheta) "
	" } "
	" # constant/default width at F-1(ALPHA) and F-1(1-ALPHA)  "
	" alpha   <- 0.1  "
	gammastr2
	lambdastr2
	hypergammastr
	hyperprior
	"} "
	];

str


% %%
% %%
model =[
	"data {"
	"  int<lower=0> N ;"
	"  int<lower=0,upper=1> y[N] ; "
	"  vector[N] x ; "
	"  real<lower=0,upper=1> alpha ; "
	fixgammastr
	fixlambdastr
	"}"
	"parameters {"
	"  real theta ;"
	"  real omega ;"
	infgammastr
	"}"
	"model {"
	" vector[N] rate;"
	" theta ~ normal(0,100);"
	" omega ~ gamma(1,0.25);"
	gammastr1
	mustr
	distr
	"}"
	];
 
model
%%

% Write the modelString to a file, using Matlab commands:
fid			= fopen(modelname,"w");
fprintf(fid,"%s\r\n",model);
fclose(fid);

%%
% keyboard
