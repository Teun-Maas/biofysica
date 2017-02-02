function samples = inferrate(x,varargin)
% MCMC = RATEFIT(X)
%
% Infer rate from promptness X
%
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
if size(x,2)>2
	warning('X should be a 1/2-column vector/matrix');
	x		= x';
end

% model

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,1000);
saveName		= keyval('saveName',varargin,'Rate-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

% graphic flags
diagFlag		= keyval('showDiag',varargin,true); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
predFlag		= keyval('showPred',varargin,true); % show posterior predictive
centroidFlag	= keyval('showCentroid',varargin,'mode'); % mode, median, mean


%% Actual regression
samples				= genMCMC(x,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic);

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
	cnt = 0;
	parameterNames	= fieldnames(samples); % get all parameter names

	n = numel(parameterNames);
	for parIdx			= 1:n
		a = samples.(parameterNames{parIdx});
		a = a(:);
		subplot(1,2,parIdx);
		cla
		plotpost(a);
		title(parameterNames{parIdx})
	end
end
% 
% %% Checks and display
% if predFlag
% 	if isnumeric(gamma),	samples.gamma	= repmat(gamma,size(samples.theta));	end
% 	if isnumeric(lambda),	samples.lambda	= repmat(lambda,size(samples.theta));	end
% 	
% 	% 	if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
% 	% 		x = x-min(x)+0.001;
% 	% 	end
% 	posteriorprediction(SOA,RT1,s,samples,fun,centroidFlag); % Posterior predictive
% 	
% 	if isnumeric(gamma),	samples		= rmfield(samples,'gamma');		end
% 	if isnumeric(lambda),	samples		= rmfield(samples,'lambda');	end
% end

%% Sub-functions
function [samples,stats] = genMCMC(x,varargin)
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

parameters		= {'muPrompt','sigmaPrompt'};

writemodel;

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.

% Specify data, as a structure
dataStruct = struct('prompt',x,'meanPrompt', mean(x),'sdPrompt',std(x),'n',numel(x));

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lambda) are determined from dependent variables/
% y data
for ii = 1:nChains
	initsStruct(ii).mu			= mean(1./x); % ~mean
	initsStruct(ii).sigma		= std(1./x); % ~mean
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



parameterNames		= fieldnames(samples) % get all parameter names


%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

function posteriorprediction(x,y,s,samples,fun,centroidFlag)
% POSTERIORPREDICTION(X,Y,B)
%
%
% keyboard
%% Parameters
theta			= samples.theta;
omega			= samples.omega;
gamma			= samples.gamma;
lambda			= samples.lambda;

% guess		= samples.guess;

%%
chainLength = length(theta);

% Graphics
us = unique(s);
ns = length(us);
sb = ceil(sqrt(ns));
if sb>5
	sb = 5;
end
figcnt	= 0;
sbcnt	= 0;
cnt		= false;

for ii = 1:ns
	
	sbcnt	= sbcnt+1;
	if cnt
		sbcnt =1;
		figcnt	= figcnt+1;
		cnt = false;
	end
	if mod(ii,25)==0
		cnt		= true;
	end
	sel		= s==us(ii);
	
	
	
	figure(200+figcnt);
	subplot(sb,sb,sbcnt)
	hold on
	% believable psychometruc curve
	cVec		= floor(linspace(1,chainLength,30));
	xWid		= max(x)-min(x);
	xComb		= linspace(min(x)-0.1*xWid,max(x)+0.1*xWid,201);
	xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
	if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}));
		selw		= xComb>=0;
		xComb	= xComb(selw);
	end
	
	% Determine centroid
	switch centroidFlag
		case 'mode' % Maximum-a-posteriori
			[f,xi]		= ksdensity(theta(:,ii));
			[~,indx]	= max(f);
			mutheta		= xi(indx);
			
			[f,xi]		= ksdensity(omega(:,ii));
			[~,indx]	= max(f);
			muomega		= xi(indx);
			
			ugamma = unique(gamma(:,ii));
			if numel(ugamma)==1
				mugamma = ugamma;
			else
				[f,xi]		= ksdensity(gamma(:,ii),'support','positive');
				[~,indx]	= max(f);
				mugamma		= xi(indx);
			end
				ulambda = unique(lambda(:,ii));
		
			if numel(ulambda)==1
				mulambda = ulambda;
			else
				[f,xi]		= ksdensity(lambda(:,ii),'support','positive');
				[~,indx]	= max(f);
				mulambda		= xi(indx);
			end
			% 				mutheta		= mode(theta(:,ii));
			% 				muomega		= mode(omega(:,ii));
			% 				mugamma		= mode(gamma(:,ii));
			% 				mulambda	= mode(lambda(:,ii));
		case 'median' % absolute error minimization
			mutheta		= median(theta(:,ii));
			muomega		= median(omega(:,ii));
			mugamma		= median(gamma(:,ii));
			mulambda	= median(lambda(:,ii));
		case 'mean' % least-squared error minimization
			mutheta		= mean(theta(:,ii));
			muomega		= mean(omega(:,ii));
			mugamma		= mean(gamma(:,ii));
			mulambda	= mean(lambda(:,ii));
	end
	
	halfway = ((1-mulambda)+mugamma)/2;
	for cIdx	= cVec
		ypred = psifun(xComb,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
		plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.7 .7 .7]);
		
		xInt	= theta(cIdx,ii);
		plot( [xInt xInt],[halfway -0.1], 'k:','Color',[.7 .7 .7]);
	end
	
	ypred = psifun(xComb,mutheta,muomega,mugamma,mulambda,0.1,'function',fun);
	
	plot(xComb,ypred,'k-','LineWidth',2,'Color','k');
	xInt		= mutheta;
	plot([xInt xInt],[halfway -0.1], 'k-','Color','k');
	horline(mugamma,'k:');
	horline(1-mulambda,'k:');
	str = {['\theta =' num2str(mutheta,'%.1f') ', \omega= ' num2str(muomega,'%.2f')],...
		['\gamma = ' num2str(round(100*mugamma)) '% , \lambda = ' num2str(round(100*mulambda)) '%']};
	text(mean(x),1.1,str,'HorizontalAlignment','center');
	
	% Data
	whos y sel
	if size(y,2)==1 % y = bernouilli rate
		[ux,~,subs]		= unique(x(sel));
		r				= accumarray(subs,y(sel),[],@sum);
		n				= accumarray(subs,ones(size(subs)),[],@sum);
	elseif size(y,2)==2 % y = [rate n]
		r = y(sel,1);
		n = y(sel,2);
		ux = x(sel);
	end
	rprop			= r./n;

	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5);
	
	%% Density
	if size(y,2)==1
		if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
			supp = 'positive';
					xi		= linspace(0.1,max(x),100);
		else
			supp = 'unbounded';
					xi		= linspace(min(x),max(x),100);
		end
		% estimate rate from Bernouilli data via ksdensity
		sel		= y==0;
		f0		= ksdensity(x(sel),xi,'support',supp);
		
		sel		= y==1;
		f1		= ksdensity(x(sel),xi,'support',supp);
		
		f = f1-f0;
		f = f-min(f);
		f = f./max(f);
		f = (1-mulambda-mugamma)*f+mugamma;
		
		plot(xi,f,'r-','LineWidth',2);
		
	end
	%% Labels
	title(us(ii))
	set(gca,'TickDir','out');
	xlabel('X');
	ylabel('Probability');
	axis square;
	box off;
	ylim([-0.1 1.2]);
	xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
	horline(halfway,'k:');
end

function writemodel
% Placeholder function to generate JAGS model

%%
str = ['# Specify the model for standardized data:\r\n',...
	'model {\r\n',...
	'\tfor ( i in 1:n ) {			\r\n',...
	'\t\t prompt[i] ~ dnorm(muPrompt,1/pow(sigmaPrompt,2))	\r\n',...
	'\t}\r\n',...
	'\t# priors\r\n',...
	'\t muPrompt ~ dnorm(meanPrompt,1/pow(sdPrompt,2))\r\n',...
	'\t sigmaPrompt ~ dunif(sdPrompt/1000,sdPrompt*1000)\r\n',...
	'}\r\n',...
	];

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);
