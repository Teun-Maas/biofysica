function b = regsat(y,x)
% B = REGJAGS(Y,X)
%
% Simple linear regression with JAGS
%
% You need to install JAGS and MATJAGS
%	http://mcmc-jags.sourceforge.net/
%	-  http://psiexp.ss.uci.edu/research/programs_data/jags/ and/or https://github.com/msteyvers/matjags
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij


%% Initialization
if nargin<2
	% Simulated height and weight data:
	[HtWtData,~,height,weight] = HtWtDataGenerator(30,5678); % from Kurschke
	x		= HtWtData(:,height);
	y		= HtWtData(:,weight);
end

%% Actual regression
b = genMCMC(x,y);

%% Pearson correlation
b.pearson = genMCMCpearson(x,y);

%% Checks and display
% posteriorprediction(x,b);

%% Sub-functions
function b = genMCMC(x,y)
% [B0,B1,SIGMA] = GENMCMC(X,Y)
%
% Generate MCMC chains
modelname = which('regress.txt');
if ~exist(modelname,'file')
	writemodel;
	modelname = fullfile(pwd, 'model.txt');
end
chain_globals;

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
[zx,mux,sdx]		= zscore(x);
[zy,muy,sdy]		= zscore(y);
nSubj				= size(x,1);

% Specify data, as a structure
dataStruct.x		= zx;
dataStruct.y		= zy;
dataStruct.Ndata	= nSubj;

%% INTIALIZE THE CHAINS.
nChains		= nChainsDefault; % Number of chains to run.
r			= corrcoef(x,y);
initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).beta0		= 0;    % because data are standardized
	initsStruct(ii).beta1		= r(2);        % because data are standardized
	initsStruct(ii).tau			= 1/(1-r(2)^2);  % because data are standardized
end

%% RUN THE CHAINS
parameters		= {'beta0' , 'beta1' , 'tau'};		% The parameter(s) to be monitored.
% adaptSteps		= 500;			% Number of steps to 'tune' the samplers.
burnInSteps		= 500;			% Number of steps to 'burn-in' the samplers.
numSavedSteps	= 5000;		% Total number of steps in chains to save.
thinSteps		= 1;			% Number of steps to 'thin' (1=keep every step).
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethodDefault,'parallel');
	doparallel		= 1; % do use parallelization
	
else
	doparallel		= 0; % do not use parallelization
end

fprintf( 'Running JAGS...\n' );
% [samples, stats, structArray] = matjags( ...
samples = matjags( ...
	dataStruct, ...                     % Observed data
	modelname, ...    % File that contains model definition
	initsStruct, ...                          % Initial values for latent variables
	'doparallel' , doparallel, ...      % Parallelization flag
	'nchains', nChains,...              % Number of MCMC chains
	'nburnin', burnInSteps,...              % Number of burnin steps
	'nsamples', nIter, ...           % Number of samples to extract
	'thin', thinSteps, ...                      % Thinning parameter
	'dic', 1, ...                       % Do the DIC?
	'monitorparams', parameters, ...     % List of latent variables to monitor
	'savejagsoutput',0, ...          % Save command line output produced by JAGS?
	'verbosity' ,1, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
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

%% Extract chain values:
z0		= [samples.beta0];
z1		= [samples.beta1];
zTau	= [samples.tau];
z0		= z0(:);
z1		= z1(:);
zTau	= zTau(:);
zSigma	= 1./zTau;

% Convert to original scale:
b1		= z1*sdy/sdx;
b0		= z0*sdy+muy-z1*sdy*mux/sdx;
sigma	= zSigma*sdy;
b.beta			= [b0 b1];
b.sigma			= sigma;


function pearson = genMCMCpearson(x,y)
%% Correlation Coefficient
chain_globals;

x		= [x y];

% Constants
[n,~]	= size(x);

%% Sampling
% MCMC Parameters
nburnin		= 500; % How Many Burn-in Samples?
nsamples	= 5000;  %How Many Recorded Samples?
nthin		= 1; % How Often is a Sample Recorded?

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n);

%% INTIALIZE THE CHAINS.
nChains		= nChainsDefault; % Number of chains to run.
initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).r		= 0;    % because data are standardized
	initsStruct(ii).mu		= zeros(1,2);        % because data are standardized
	initsStruct(ii).lambda = ones(1,2);  % because data are standardized
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
	'verbosity' ,1, ...
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

pearson = struct('r',r','mu',mu,'sigma',sigma);


function posteriorprediction(x,b)
% POSTERIORPREDICTION(X,Y,B)
%
%

b0		= b.beta(:,1);
b1		= b.beta(:,2);
sigma	= b.sigma;

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


function writemodel
% Placeholder function to generate JAGS model
% y = 90*2*(logistic(x/90*4)-0.5);

str = ['model {\r\n',...
	'\tfor( i in 1 : Ndata ) {\r\n',...
	'\t\ty[i] ~ dnorm( mu[i] , tau )\r\n',...
	'\t\tmu[i] <- 90*2*(ilogit( (beta0 + beta1 * x[i])/90*4)-0.5) )\r\n',...
	'\t}\r\n',...
	'\tbeta0 ~ dnorm( 0 , 1.0E-12 )\r\n',...
	'\tbeta1 ~ dnorm( 0 , 1.0E-12 )\r\n',...
	'\ttau ~ dgamma( 0.001 , 0.001 )\r\n',...
	'}\r\n',...
	];

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);

