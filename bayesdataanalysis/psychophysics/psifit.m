function samples = psifit(x,y,s,varargin)
% MCMC = PSYCHCURVE(X,Y,S)
%
% Psychometric curve
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



%% Initialization
if isempty(s)
	s = ones(1,size(y,2));
end
% check for input orientation
% Y needs to be either:
% - a nx1 column vector with 0s and 1s
% - a nx2 column matrix with rate in colum 1 and #trials in column 2
if size(y,2)>2
	warning('Y should be a 1/2-column vector/matrix');
	y = y';
	x = x';
	s = s';
end
if size(y,1)~=size(x,1)
	warning('Y and X should have same orientations');
	x = x';
end
if size(y,1)~=size(s,1)
	warning('Y and S should have same orientations');
	s = s';
end
if size(y,1)~=size(x,1)
	error('X and Y should have same number of rows');
end
%% S should contain consecutive values, no missing values
[~,~,s] = unique(s);

% model
if size(y,2)==2
	respDist = 'binomial';
elseif size(y,2)==1
	respDist = 'bernouilli';
end


fun				= keyval('function',varargin,@logisticfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
fig				= keyval('figure',varargin,200); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun

gamma			= keyval('gamma',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 gamma for all subjects/conditions
lambda			= keyval('lambda',varargin,'infer');% 0.001 = low lapse rate, 'infer' = estimate from data, 'single' = infer 1 lambda for all subjects/conditions

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,1000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

% graphic flags
diagFlag		= keyval('showDiag',varargin,false); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
predFlag		= keyval('showPred',varargin,true); % show posterior predictive
centroidFlag	= keyval('showCentroid',varargin,'mean'); % mode, median, mean

%% Actual regression
samples				= genMCMC(x,y,s,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,'function',fun,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic,'gamma',gamma,'lambda',lambda,'respDist',respDist);

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
	parameterNames	= {'mutheta','muomega','mugamma','mulambda'}; % get all parameter names
	n = numel(parameterNames);
	for parIdx1			= 1:n
		for parIdx2			= 1:n
			cnt = cnt+1;
			a = samples.(parameterNames{parIdx1});
			b = samples.(parameterNames{parIdx2});
			
			subplot(n,n,cnt)
			plot(a,b,'.');
			box off
			axis square
			xlabel(parameterNames{parIdx1});
			ylabel(parameterNames{parIdx2});
		end
	end
	
	for parIdx			= 1:n
		a = samples.(parameterNames{parIdx});
		a = a(:);
		subplot(n,n,n*(parIdx-1)+parIdx);
		cla
		plotpost(a);
		
	end
end

%% Checks and display

if predFlag
	if isnumeric(gamma),	samples.gamma	= repmat(gamma,size(samples.theta));	end
	if isnumeric(lambda),	samples.lambda	= repmat(lambda,size(samples.theta));	end
	figure(fig);

	% 	if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	% 		x = x-min(x)+0.001;
	% 	end
	posteriorprediction(x,y,s,samples,fun,centroidFlag,fig); % Posterior predictive
	
	if isnumeric(gamma),	samples		= rmfield(samples,'gamma');		end
	if isnumeric(lambda),	samples		= rmfield(samples,'lambda');	end
end

%% Sub-functions
function [samples,stats] = genMCMC(x,y,s,varargin)
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
fun				= keyval('function',varargin);
runjagsMethod	= keyval('runjagsMethod',varargin);
dic				= keyval('dic',varargin);
gamma			= keyval('gamma',varargin);
lambda			= keyval('lambda',varargin);
modelname		= fullfile(pwd, 'model.txt');
respDist		= keyval('respDist',varargin);

%% Write the model
ngroups = numel(unique(s));

% first check parameters to be monitored
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
writemodel(gamma,lambda,ngroups,fun,respDist);

% Re-center data at mean, to reduce autocorrelation in MCMC sampling.
% Standardize (divide by SD) to make initialization easier.
[zx,mux,sdx]		= zscore(x);
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	% Weibull functions are defined for x>0
	mux = min(x);
	% 	mux = 0;
	sdx = 1;
	zx	= (x-mux)/sdx;
end
Ntotal				= size(x,1);
Nsubj				= max(s);
% Specify data, as a structure
dataStruct = struct('x',zx,'y',y,'s',s,'Ntotal',Ntotal,'Nsubj',Nsubj);

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lambda) are determined from dependent variables/
% y data
for ii = 1:nChains
	initsStruct(ii).theta			= zeros(Nsubj,1); % ~mean
	initsStruct(ii).omega		= ones(Nsubj,1);  % ~std
	if ~isnumeric(gamma)
		initsStruct(ii).gamma		= repmat(0.01,Nsubj,1); % no guessing
	end
	if ~isnumeric(lambda)
		initsStruct(ii).lambda		= repmat(0.01,Nsubj,1); % no lapses
	end
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



parameterNames		= fieldnames(samples); % get all parameter names
for ii					= 1:numel(parameterNames)
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



%% Save samples
if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

function posteriorprediction(x,y,s,samples,fun,centroidFlag,fig)
% POSTERIORPREDICTION(X,Y,B)
%
%

%% Parameters
theta			= samples.theta;
omega			= samples.omega;
gamma			= samples.gamma;
lambda			= samples.lambda;

% guess		= samples.guess;

%% Determine centroid
if strcmp(centroidFlag,'mode')
	centroid = ['bf_' centroidFlag]; % for random samples, BF_MODE works better than MODE
else
	centroid = centroidFlag;
end
parameterNames	= fieldnames(samples); % get all parameter names
for parIdx			= 1:numel(parameterNames)
	str = [parameterNames{parIdx} ' = samples.(parameterNames{parIdx});'];
	eval(str);
	str = ['c' parameterNames{parIdx} ' = ' centroid '(samples.(parameterNames{parIdx}));'];
	eval(str);
end

%%
chainLength = length(theta);

%% believable psychometric curve
cVec		= floor(linspace(1,chainLength,30));
xWid		= max(x)-min(x);
xComb		= linspace(min(x)-0.1*xWid,max(x)+0.1*xWid,201);
xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	selw		= xComb>=0;
	xComb		= xComb(selw);
end

%% some shenangins to plot a max of 25 subplots per figure
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
	% for every subject/group
	sbcnt		= sbcnt+1; % new subplot count
	if cnt
		sbcnt	= 1;
		figcnt	= figcnt+1;
		cnt		= false;
	end
	if mod(ii,25)==0 % when 25 subplots have been reached
		cnt		= true; % new figure (see above)
	end
	sel			= s==us(ii);
	
	figure(fig+figcnt);
	subplot(sb,sb,sbcnt)
	hold on
	
	%% Posterior predictive
	% Data
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
	
	%
% 	for cIdx	= cVec
% 		ypred	= psifun(ux,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
% 		whos n ypred
% 		keyboard
% 		R = binornd(n,ypred);
% 	end
	
	%% Posterior credible psychometric curves
	halfway = ((1-clambda)+cgamma)/2;
	for cIdx	= cVec
		ypred	= psifun(xComb,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
		xInt	= theta(cIdx,ii);

		plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.7 .7 .7]);
		plot( [xInt xInt],[halfway(ii) -0.1], 'k:','Color',[.7 .7 .7]);
	end

	
	%% Max Posterior Curve
	ypred = psifun(xComb,ctheta(ii),comega(ii),cgamma(ii),clambda(ii),0.1,'function',fun);
	
	
	%% Graphics
	plot(xComb,ypred,'k-','LineWidth',2,'Color','k');
	xInt		= ctheta(ii);
	plot([xInt xInt],[halfway(ii) -0.1], 'k-','Color','k');
	horline(cgamma(ii),'k:');
	horline(1-clambda(ii),'k:');
	str = {['\theta =' num2str(ctheta(ii),'%.1f') ', \omega= ' num2str(comega(ii),'%.2f')],...
		['\gamma = ' num2str(round(100*cgamma(ii))) '% , \lambda = ' num2str(round(100*clambda(ii))) '%']};
	text(mean(x),1.1,str,'HorizontalAlignment','center');
	
    
	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5); % data
	
% 	%% Density
% 	if size(y,2)==1
% 		if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
% 			supp = 'positive';
% 			xi		= linspace(0.1,max(x),100);
% 		else
% 			supp = 'unbounded';
% 			xi		= linspace(min(x),max(x),100);
% 		end
% 		% estimate rate from Bernouilli data via ksdensity
% 		sel		= y==0;
% 		f0		= ksdensity(x(sel),xi,'support',supp);
% 		
% 		sel		= y==1;
% 		f1		= ksdensity(x(sel),xi,'support',supp);
% 		
% 		f = f1-f0;
% 		f = f-min(f);
% 		f = f./max(f);
% 		f = (1-clambda(ii)-cgamma(ii))*f+cgamma(ii);
% 		
% 		plot(xi,f,'r-','LineWidth',2);
% 		
% 	end
	%% Labels
	title(us(ii))
	set(gca,'TickDir','out');
	xlabel('X');
	ylabel('Probability');
	axis square;
	box off;
	ylim([-0.1 1.2]);
	xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);
	horline(halfway(ii),'k:');
end

function writemodel(gamma,lambda,ngroups,fun,respDist)
% Placeholder function to generate JAGS model

%% Psychometric function
fun = func2str(fun);
switch fun
	case 'logisticfun'
		fun = 'ilogit';
		par = '(2*log(1/alpha-1))/omega[s[i]] * (x[i]-theta[s[i]])';
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
		fun = 'icloglog';
		par = '( (2*omega[s[i]]*theta[s[i]])/log(2) * (log(x[i])-log(theta[s[i]])) + log(log(2)) )';
	case 'revweibullfun'
		fun = '1-icloglog';
		par = '( (-2*omega[s[i]]*theta[s[i]])/log(2) * (log(x[i])-log(theta[s[i]])) + log(log(2)) )';
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
hyperprior = ['\t# threshold/theta hyperprior \r\n',...
	'\t mutheta\t\t\t~ dnorm(0,1/(10)^2)\r\n',...
	'\t sigmatheta\t\t~ dunif(1.0E-3,1000)\r\n',...
	'\t# width/omega hyperprior \r\n',...
	'\t somega\t\t\t<- pow(muomega,2)/pow(sigmaomega,2) # shape \r\n',...
	'\t romega\t\t\t<- muomega/pow(sigmaomega,2) # rate \r\n',...
	'\t muomega\t\t\t\t~ dgamma(1,0.25)\r\n',...
	'\t sigmaomega\t\t\t\t~ dgamma(1,0.5) \r\n'];
else
hyperprior = ['\t# threshold/theta hyperprior \r\n',...
	'\t mutheta\t\t\t<- 0\r\n',...
	'\t sigmatheta\t\t<- 1/(10)^2\r\n',...
	'\t# width/omega hyperprior \r\n',...
	'\t somega\t\t\t<- 1.0 # shape \r\n',...
	'\t romega\t\t\t<- 0.25 # rate \r\n'];
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
	'\t\t omega[j]\t~ dgamma(somega,romega)  \r\n',...
	'\t\t theta[j]\t~ dnorm(mutheta,sigmatheta)\r\n',...
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
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);
