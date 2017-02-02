function samples = twinrawfit(t,RT,RT1,RT2,varargin)
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
burnInSteps		= keyval('burnInSteps',varargin,1000);
saveName		= keyval('saveName',varargin,'TWIN-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

% graphic flags
diagFlag		= keyval('showDiag',varargin,false); % show MCMC diagnostics
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

parameters		= {'omega','mu','delta','kappa','gamma','sigma','PI','PW','lambda1mu','lambda2mu'};

writemodel;

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
nt					= numel(t,1);

% Specify data, as a structure
dataStruct = struct('t',t,'rt',rt,'nt',nt,...
	'prompt1',1./rt1,'nrt1',numel(rt1),'meanprompt1',mean(1./rt1),'sdprompt1',std(1./rt1),...
	'prompt2',1./rt2,'nrt2',numel(rt2),'meanprompt2',mean(1./rt2),'sdprompt2',std(1./rt2));

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
	initsStruct(ii).sigma		= 10; % ~mean
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
% 
% d = 142;
% m = 250;
% o = 342;
% lv = 1/220;
% la = 1/45;
% ga = 300;
% k = 20;
% 
% ypred		= twinfun(xComb,lv,la,m,o,d,ga,k);
% plot(xComb, ypred,'k-','LineWidth',1.5,'Color','k');
% % RT = twinfun(T,lv,la,m,o,d,ga,k);
% 
% %
% 	RT = twinfun(T,lambda1,lambda2,omega,delta,gamma,kappa)
for cIdx	= cVec
	ypred = twinfun(xComb,lambda1mu(cIdx),lambda2mu(cIdx),mu(cIdx),omega(cIdx),delta(cIdx),gamma(cIdx),kappa(cIdx));
	plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.7 .7 .7]);
end


% 	ypred = psifun(xComb,mutheta,muomega,mugamma,mulambda,0.1,'function',fun);
%
% 	plot(xComb,ypred,'k-','LineWidth',2,'Color','k');
% 	xInt		= mutheta;
% 	plot([xInt xInt],[halfway -0.1], 'k-','Color','k');
% 	horline(mugamma,'k:');
% 	horline(1-mulambda,'k:');
% 	str = {['\theta =' num2str(mutheta,'%.1f') ', \omega= ' num2str(muomega,'%.2f')],...
% 		['\gamma = ' num2str(round(100*mugamma)) '% , \lambda = ' num2str(round(100*mulambda)) '%']};
% 	text(mean(t),1.1,str,'HorizontalAlignment','center');
%
% 	% Data
% 	whos y sel
% 	if size(rt,2)==1 % y = bernouilli rate
% 		[ux,~,subs]		= unique(t(sel));
% 		r				= accumarray(subs,rt(sel),[],@sum);
% 		n				= accumarray(subs,ones(size(subs)),[],@sum);
% 	elseif size(rt,2)==2 % y = [rate n]
% 		r = rt(sel,1);
% 		n = rt(sel,2);
% 		ux = t(sel);
% 	end
% 	rprop			= r./n;
%
% 	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5);
%
% 	%% Density
% 	if size(rt,2)==1
% 		if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
% 			supp = 'positive';
% 					xi		= linspace(0.1,max(t),100);
% 		else
% 			supp = 'unbounded';
% 					xi		= linspace(min(t),max(t),100);
% 		end
% 		% estimate rate from Bernouilli data via ksdensity
% 		sel		= rt==0;
% 		f0		= ksdensity(t(sel),xi,'support',supp);
%
% 		sel		= rt==1;
% 		f1		= ksdensity(t(sel),xi,'support',supp);
%
% 		f = f1-f0;
% 		f = f-min(f);
% 		f = f./max(f);
% 		f = (1-mulambda-mugamma)*f+mugamma;
%
% 		plot(xi,f,'r-','LineWidth',2);
%
% 	end
% 	%% Labels
% 	title(us(ii))
% 	set(gca,'TickDir','out');
% 	xlabel('X');
% 	ylabel('Probability');
% 	axis square;
% 	box off;
% 	ylim([-0.1 1.2]);
% 	xlim([min(t)-0.1*xWid,max(t)+0.1*xWid]);
% 	horline(halfway,'k:');
% end

function writemodel
% Placeholder function to generate JAGS model

%%
str = ['# Specify the model for standardized data:\r\n',...
	'model {\r\n',...
	'\t# Bimodal\r\n',...
	'\tfor ( i in 1:nt ) {			\r\n',...
	'\t\t rt[i] ~ dnorm(rtmu[i],1/pow(sigma,2))	\r\n',...
	'\t\t rtmu[i] <- 1/lambda1mu+mu-delta*PI[i]-kappa*PW[i]		\r\n',...
	'\t\t# for every stimulus onset asynchrony t\r\n',...
	'\t\t# Probability of integration P(I)\r\n',...
	'\t\t PI[i] <- ifelse(x1[i],a[i],ifelse(x2[i],b[i],c[i]))		\r\n',...
	'\t\t x1[i] <- t[i]<(t[i]+omega) && (t[i]+omega)<0	\r\n',...
	'\t\t x2[i] <- t[i]<=0 && 0<=(t[i]+omega)		\r\n',...
	'\t\t x3[i] <- 0<t[i] && t[i]<(t[i]+omega)		\r\n',...
	'\t\t a[i] <- lambda1mu/(lambda1mu+lambda2mu)*exp(lambda2mu*t[i])*(-1+exp(lambda2mu*omega))				\r\n',...
	'\t\t b[i] <- 1/(lambda1mu+lambda2mu)*(lambda2mu*(1-exp(-lambda1mu*(omega+t[i])))+lambda1mu*(1-exp(lambda2mu*t[i]))) \r\n',...
	'\t\t c[i] <- lambda2mu/(lambda1mu+lambda2mu)*(exp(-lambda1mu*t[i])-exp(-lambda1mu*(omega+t[i])))			\r\n',...
	'\t\t# Probability of warning P(W)\r\n',...
	'\t\t PW[i] <- ifelse(x4[i],d[i],e[i])		\r\n',...
	'\t\t x4[i] <- (t[i]+gamma)<0	\r\n',...
	'\t\t d[i] <- 1-(lambda1mu/(lambda1mu+lambda2mu)*exp(lambda2mu*(t[i]+gamma)))			\r\n',...
	'\t\t e[i] <- (lambda2mu/(lambda1mu+lambda2mu)*exp(-lambda1mu*(t[i]+gamma)))	\r\n',...
	'\t}\r\n',...
	'\t# Bimodal priors\r\n',...
	'\tomega ~ dgamma(1,1)\r\n',...
	'\tdelta ~ dgamma(1,1)\r\n',...
	'\tkappa ~ dgamma(1,1)\r\n',...
	'\tgamma ~ dgamma(1,1)\r\n',...
	'\t# Unimodal\r\n',...
	'\t\t lambda1mu <- 1/(1/prompt1mu-mu)			\r\n',...
	'\t\t lambda2mu <- 1/(1/prompt2mu-mu)			\r\n',...
	'\tfor ( j1 in 1:nrt1 ) {			\r\n',...
	'\t\t prompt1[j1] ~ dnorm(prompt1mu,1/pow(prompt1sigma,2))	\r\n',...
	'\t}\r\n',...
	'\tfor ( j2 in 1:nrt2 ) {			\r\n',...
	'\t\t prompt2[j2] ~ dnorm(prompt2mu,1/pow(prompt2sigma,2))	\r\n',...
	'\t}\r\n',...
	'\t# Unimodal priors\r\n',...
	'\tmu ~ dgamma(1,1)\r\n',...
	'\tsigma ~ dgamma(1,1)\r\n',...
	'\tprompt1mu ~ dnorm(meanprompt1,1/pow(sdprompt1,2))\r\n',...
	'\tprompt2mu ~ dnorm(meanprompt2,1/pow(sdprompt2,2))\r\n',...
	'\tprompt1sigma ~ dgamma(1,1)\r\n',...
	'\tprompt2sigma ~ dgamma(1,1)\r\n',...
	'}\r\n',...
	];
% 	'\tl1 ~ dgamma(0.001,0.001)\r\n',...
% 	'\tl2 ~ dgamma(0.001,0.001)\r\n',...

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);
