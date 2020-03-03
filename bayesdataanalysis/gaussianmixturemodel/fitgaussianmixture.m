function [samples,stats] = fitgaussianmixture(data,Nclust,varargin)
% SAMPLES = FITGAUSSIANMIXTURE(DATA,NCLUST)
%
% Fit a mixture of NCLUST number of gaussian distributions on DATA via
% Bayesian MCMC methods, yielding posterior SAMPLES on:
% - SAMPLES.mu		= mean of cluster
% - SAMPLES.sigma	= standard deviation of cluster 
% - SAMPLES.p		= proportion of data samples in cluster
%
% SAMPLES also contains a posterior distribution of samples indicatinng to
% which cluster each datum value was assigned to, SAMPLES.clust.
%
% SAMPLES = FITGAUSSIANMIXTURE(DATA,NCLUST,'KEY',VALUE)
%
% - key: 'CLUST', value: a NaN(N,1) vector for which at least NCLUST number
% of data points have each been assigned to a unique cluster.
%
% For the clustering, you need to assign at least one data point to each
% cluster. By default, the data point closest to 0 will be set to belong to
% cluster 1, the largest data point will be set to belong to cluster 2
% (this corresponds to one cluster belonging to a variable that is
% independent of any experimental manipulation, while the other is
% dependent).
%
% Other key-value inputs are the default sampling parameters, such as:
% - numSavedSteps (default 1000)
% - burnInSteps (default 500)
%
% See also: MATJAGS, CHAIN_GLOBALS, PLOTPOST
%
% See also: http://doingbayesiandataanalysis.blogspot.com/2012/06/mixture-of-normal-distributions.html


%% Initialization
if nargin<1
% Generate random data from known parameter values:
	
	% 	trueM1		= 100;
	% 	N1			= 500;
	% 	trueM2		= 145; % 145 for first example below; 130 for second example
	% 	N2			= 500;
	% 	trueSD		= 15;
	
	trueM1		= 0;
	trueM2		= 0;
	N1			= 200;
	N2			= 300;
	trueSD1		= 1;
	trueSD2		= 10;
	
	y1			= randn(N1,1);
	y1			= (y1-mean(y1))/std(y1) * trueSD1 + trueM1;
	y2			= randn(N2,1);
	y2			= (y2-mean(y2))/std(y2) * trueSD2 + trueM2;
	data		= [y1; y2];
end
N				= length(data);

if nargin<2
	Nclust		= 2;
end

%% Initial clustering
% Must have at least one data point with fixed assignment
% to each cluster, otherwise some clusters will end up empty:
clust			= NaN(N,1);
[~,mn]			= min(abs(data));
[~,mx]			= max(data);
clust(mn)		= 1; % smallest value assigned to cluster 1
clust(mx)		= 2; % highest value assigned to cluster 2
clust			= keyval('clust',varargin,clust);

%% matJAGS parameters
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,1000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,500);
saveName		= keyval('saveName',varargin,'Hier-GMM-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end

%% Other KEY-VALUE inputs
if nargout<1
	graph			= keyval('showPlot',varargin,true);
else
	graph			= keyval('showPlot',varargin,false);
end

%% Data
dataStruct		= struct('y',data,....
	'N',N,'Nclust',Nclust,'clust',clust,'onesRepNclust',ones(Nclust,1)); % data structure
modelname		= which('fitgmm.txt');
parameters		= {'muOfClust','tauOfClust','pClust','clust'};

%% Initial guess
for ii = 1:nChains
	initsStruct(ii).muOfClust			= [0 mean(data)]; % standardized mean
	initsStruct(ii).tauOfClust			= [1 1./std(data).^2]; % standardized mean
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


%% Vectorizing chain samples
samples = extractchain(samples);

%% Rename fieldss
samples			= replacexy(samples,{'mu','sigma','p'},{'muOfClust','tauOfClust','pClust'});
samples.sigma	= 1./sqrt(samples.sigma);
%% Graphics
if graph
	pC1 = samples.p(:,1);
	pC2 = samples.p(:,2);
	mC2 = samples.mu(:,2);
	mC1 = samples.mu(:,1);
	sC1 = samples.sigma(:,1);
	sC2 = samples.sigma(:,2);
	
	c = mean(samples.clust)-1;
	
	figure(1)
	clf
	subplot(231)
	
	[f,x] = ksdensity(mC1);
	hold on
	
	plot(x,f,'LineWidth',2);
	hold on
	[f,x] = ksdensity(mC2);
	plot(x,f,'LineWidth',2);
	nicegraph;
	xlabel('mean');
	if exist('trueM1','var')
	verline([trueM1 trueM2]);
	end
	subplot(232)
	[f,x] = ksdensity(sC1);
	plot(x,f,'-','LineWidth',2);
	
	hold on
	[f,x] = ksdensity(sC2);
	plot(x,f,'-','LineWidth',2);
	
	nicegraph;
	xlabel('standard deviation');
	if exist('trueM1','var')
		verline([trueSD1 trueSD2]);
	end
	
	subplot(233)
	[f,x] = ksdensity(pC1);
	plot(x,f,'LineWidth',2);
	hold on
	[f,x] = ksdensity(pC2);
	plot(x,f,'LineWidth',2);
	nicegraph;
	xlabel('proportion');
	verline([N1 N2]./(N1+N2));
	
	subplot(234)
	plot(mC1,sC1,'.');
	hold on
	plot(mC2,sC2,'.');
	nicegraph;
	xlabel('mean');
	ylabel('standard deviation');
	
	subplot(236)
	[~,idx] = sort(data);
	plot(data(idx),c(idx),'k-','MarkerFaceColor','w','LineWidth',2);
	nicegraph;
	ylabel('probability cluster 2');
	xlabel('datum value');
	
	subplot(235)
	plotpost(data(:));
end

%% Check whether output is asked for
if nargout<1
	clear samples stats
end