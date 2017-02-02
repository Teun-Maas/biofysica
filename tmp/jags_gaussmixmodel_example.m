function jags_gaussmixmodel_example

%% Initialization
close all

%% The data specification:
 [y,N] = getdata;

%% JAGS model
[samples,nChains] = jagsmodel(y,N);

%% Check samples
samples = clustering(samples,nChains);
checkchaincluster(samples);
samples = extractchain(samples);

%% Plot relevant parameters
C = NaN(N,1);
for ii = 1:N
	C(ii) = mean(samples.C(:,ii)==1);
end

figure(3)
clf
subplot(221)
plot(y)

subplot(222)
hist(y,0:5:200)
xlim([0 200])

subplot(223)
plot(C)
ylim([-0.1 1.1])

subplot(224)
plotpost(samples.mu(:,1))

plotpost(samples.mu(:,2))
xlim([0 200])


%%

function [y,N] = getdata
% Generate random data from known parameter values:
trueM1	= 100;
N1		= 400;
trueM2	= 150; % 145 for first example below; 130 for second example
N2		= 200;
trueSD	= 5;
y1		= randn(N1,1);
y1 = (y1-mean(y1))/std(y1) * trueSD + trueM1;
y2 = randn(N2,1);
y2 = (y2-mean(y2))/std(y2) * trueSD + trueM2;
y = [y1; y2];
N = length(y);

function [samples,nChains] = jagsmodel(y,N)
chain_globals;
dataStruct = struct('d',y,...
	'nobs',N);
runjagsMethod = runjagsMethodDefault;
if strcmp(runjagsMethod,'parallel');
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
modelname	= which('mixmodel.txt');
nChains		= 3;
burnInSteps = 1000;
nIter		= 1000;
thinSteps	= 1;
dic			= 0;
parameters	= {'C','pi','mu','tau'};
initsStruct = struct([]);
for ii		= 1:nChains
	initsStruct(ii).tau		= [10 10];  % ~std
end
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
	'dic',dic, ...                       % Do the DIC?
	'monitorparams', parameters, ...     % List of latent variables to monitor
	'savejagsoutput',0, ...          % Save command line output produced by JAGS?
	'verbosity',0, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
	'cleanup',1);                    % clean up of temporary files?

function samples = clustering(samples,nChains)

%% Get chain clusters
mu1 = mean(squeeze(samples.mu(:,:,1)),2);
mu2 = mean(squeeze(samples.mu(:,:,2)),2);
idx = [1 1 1; 1 1 0; 1 0 1; 0 1 1];
sd = NaN(4,1);
for ii = 1:4
	sel = logical(idx(ii,:));
	sd(ii) = std([mu1(sel); mu2(~sel)]);
end
[~,sdidx]	= min(sd);
idx			= idx(sdidx,:);

%% Rearrange
parameterNames	= fieldnames(samples); % get all parameter names
for parIdx			= 1:numel(parameterNames)
	for ii = 1:nChains
		if idx(ii)
					n = size(samples.(parameterNames{parIdx}),3);
					if n==2

			a = samples.(parameterNames{parIdx})(ii,:,1);
			b = samples.(parameterNames{parIdx})(ii,:,2);
			samples.(parameterNames{parIdx})(ii,:,1) = b;
			samples.(parameterNames{parIdx})(ii,:,2) = a;
					else
									c = samples.(parameterNames{parIdx})(ii,:,:);
									sel = c==1;
									c(sel)=2;
									c(~sel) = 1;
 samples.(parameterNames{parIdx})(ii,:,:) =  c;
					end
			
		end
	end
end

function checkchaincluster(samples)

%% Check clustering
figure(2)
clf
subplot(221)
mu = squeeze(samples.mu(:,:,1));
plot(mu');

subplot(222)
mu = squeeze(samples.mu(:,:,2));
plot(mu');

subplot(223)
mu = squeeze(samples.pi(:,:,1));
plot(mu');

subplot(224)
mu = squeeze(samples.pi(:,:,2));
plot(mu');
