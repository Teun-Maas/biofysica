function samples = exppsifit(x,y,s,varargin)
% MCMC = PSIFIT(X,Y,S)


%% Initialization

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
if size(y,1)~=size(x,1)
	error('X and Y should have same number of rows');
end
if isempty(s)
	s = ones(size(y,1),1);
end
%% S should contain consecutive values, no missing values
[~,~,s] = unique(s);

% model
if size(y,2)==2
	respDist = 'binomial';
elseif size(y,2)==1
	respDist = 'bernouilli';
end


fig				= keyval('figure',varargin,200); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
alpha			= keyval('alpha',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 alpha for all subjects/conditions

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
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

	samples = jags_expfit(x,y,s,...
		'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
		'saveName',saveName,'nChains',nChains,...
		'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
		'dic',dic,'alpha',alpha,'respDist',respDist);


%% Thin
if thinSteps>1
	samples = thinchain(samples,thinSteps);
end

%% MCMC diagnostics
if diagFlag
	diagmcmc(samples)
end

%% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

%% Posterior estimates
if postFlag
	figure;
	cnt = 0;
	parameterNames	= {'mulambda','mualpha','mulambda'}; % get all parameter names
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

	figure(fig);

	% 	if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}))
	% 		x = x-min(x)+0.001;
	% 	end
	posteriorprediction(x,y,s,samples,centroidFlag,fig); % Posterior predictive
	
	
end


function posteriorprediction(x,y,s,samples,centroidFlag,fig)
% POSTERIORPREDICTION(X,Y,B)
%
%

%% Parameters
alpha			= samples.alpha;
lambda			= samples.lambda;

% guess		= samples.guess;

%% Determine centroid
if strcmp(centroidFlag,'mode')
	centroid = ['bf_' centroidFlag]; % for random samples, BF_MODE works better than MODE
else
	centroid = centroidFlag;
end
parameterNames	= fieldnames(samples) % get all parameter names
for parIdx			= 1:numel(parameterNames)
	str = [parameterNames{parIdx} ' = samples.(parameterNames{parIdx});'];
	eval(str);
	str = ['c' parameterNames{parIdx} ' = ' centroid '(samples.(parameterNames{parIdx}));'];
	eval(str);
end

%%
chainLength = length(lambda)

%% believable psychometric curve
cVec		= floor(linspace(1,chainLength,30));
xWid		= max(x)-min(x);
xComb		= linspace(min(x),max(x),201);
xlim([min(x)-0.1*xWid,max(x)+0.1*xWid]);


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
	
	% Posterior predictive
	% Data
	if size(y,2)==1 % y = bernouilli rate
		
		[r,ux]	= assemble(y(sel),x(sel),'fun',@sum);
		n		= assemble(ones(size(y(sel))),x(sel),'fun',@sum);
	elseif size(y,2)==2 % y = [rate n]
		[r,ux]	= assemble(y(sel,1),x(sel),'fun',@sum);
		n		= assemble(y(sel,2),x(sel),'fun',@sum);		
	end
	rprop			= r./n;
	[lb, ub]		= binomialci(r, n, 0.05);	
	
	
% Posterior credible psychometric curves
	halfway = ((1-clambda)+calpha)/2;
	for cIdx	= cVec
		ypred = alpha(cIdx,ii) * exp( -1/(lambda(cIdx,ii))^2 * xComb.^2 );
		xInt	= lambda(cIdx,ii);

		plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.7 .7 .7]);
		plot( [xInt xInt],[halfway(ii) -0.1], 'k:','Color',[.7 .7 .7]);
	end

	
	% Max Posterior Curve
% 	ypred = psifun(xComb,clambda(ii),calpha(ii),clambda(ii),'function',fun);
% 	ypred = psifun(xComb,clambda(ii),calpha(ii),clambda(ii),'function',fun);
	
	ypred = calpha(ii) * exp( -1/clambda(ii)^2 * xComb.^2 );
	
	% Graphics
	plot(xComb,ypred,'k-','LineWidth',2,'Color','k');
	xInt		= clambda(ii);
	plot([xInt xInt],[halfway(ii) -0.1], 'k-','Color','k');
	horline(calpha(ii),'k:');
	horline(1-clambda(ii),'k:');
	str = {['\lambda =' num2str(clambda(ii),'%.1f') ],...
		['\alpha = ' num2str(calpha(ii),'%0.2f') ', \lambda = ' num2str(clambda(ii),'%0.1f')]};
	text(mean(x),1.1,str,'HorizontalAlignment','center');
	
    verline(clambda(ii));
	horline(calpha(ii)*exp(-1));
% 	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5); % data
	errorbar(ux,rprop,rprop-lb,ub-rprop,'ks','MarkerFaceColor','w','MarkerSize',10); % data
	
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
% 		f = (1-clambda(ii)-calpha(ii))*f+calpha(ii);
% 		
% 		plot(xi,f,'r-','LineWidth',2);
% 		
% 	end
	% Labels
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




function [samples,stats] = jags_expfit(x,y,s,varargin)
% SAMPLES = JAGS_PSIFIT(X,Y)
%
% See also PSIFIT

%% Psychometric function
fun				= keyval('function',varargin,@expfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
alpha			= keyval('alpha',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 alpha for all subjects/conditions
% modelname		= fullfile(pwd, 'jags_exp_model.txt');
modelname = '/Users/marcw/Dropbox/Manuscript/Luuk van de Rijt/#4 SPIN AV CI/matlab/exp_model.txt';
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
% writemodel(alpha,Ngroups,fun,respDist,modelname);

%% Data
parameters		= {'lambda','alpha','mulambda','sigmalambda','mugamma','kappagamma'};
dataStruct		= struct('x',x,'y',y,'s',s,'Ntotal',Ntotal,'Nsubj',Ngroups); % data structure
initsStruct		= initialize(Ngroups,nChains); % initial parameter values

%% parallel?
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
% 	doparallel		= 0; % do not use parallelization

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


function initsStruct	= initialize(Nsubj,nChains)
initsStruct = struct([]);
% because covariates/explanatory variables/x data are standardized,
% means(x) can be set to zero, and stds(x) to 1.
% Probabilities (guess and lapse rate) are determined from dependent variables/
% y data
for ii = 1:nChains
	initsStruct(ii).lambda			=repmat(20,Nsubj,1); % standardized mean
		initsStruct(ii).alpha		=repmat(0.99,Nsubj,1); % no guessing
end

function samples		= unzscore(samples,mux,sdx)
parameterNames		= fieldnames(samples); % get all parameter names
for ii				= 1:numel(parameterNames)
	z			= samples.(parameterNames{ii});
	% Convert to original scale:
	switch parameterNames{ii}
		case 'lambda'
			samples.(parameterNames{ii})		= z*sdx+mux;
		case 'mulambda'
			samples.(parameterNames{ii})		= z*sdx+mux;

	end
end

function writemodel(alpha,ngroups,fun,respDist,modelname)
% Placeholder function to generate JAGS model


%% Psychometric function

		fun = 'exp';
		par = '(-x[i]*lambda[s[i]])';

%% Guess rate alpha
if isnumeric(alpha)
	alphastrmu	= 'alpha';
	alphastr1	= 	[];
	alphastr2	= 	['\t alpha\t\t\t<- ' num2str(alpha) '  \r\n'];
elseif strcmp(alpha,'infer')
	alphastrmu = 'alpha[s[i]]';
	alphastr1	=	'\t\t alpha[j]\t~ dbeta(aalpha,balpha)  \r\n';
	alphastr2	= [	'\t# guess rate / alpha hyperprior \r\n',...
		'\t aalpha\t\t\t<- mualpha*kappaalpha \r\n',...
		'\t balpha\t\t\t<- (1.0-mualpha) * kappaalpha  \r\n',...
		'\t mualpha\t\t~ dbeta(2.0,2.0)  \r\n',...
		'\t kappaalpha\t\t~ dalpha(Skappa,Rkappa)  \r\n'];
	if ngroups==1
			alphastr2	= [	'\t# guess rate / alpha hyperprior \r\n',...
		'\t aalpha\t\t\t<- 1.0 \r\n',...
		'\t balpha\t\t\t<- 1.0  \r\n'];
	end
end

%% Lapse rate lambda


if isempty(alphastr1)
	hyperalphastr = [];
else
	hyperalphastr = [	'\t# constants for all dalpha hyperprior \r\n',...
		'\t Skappa\t\t\t<- pow(10,2)/pow(10,2) # constant \r\n',...
		'\t Rkappa\t\t\t<- 10/pow(10,2) # constant \r\n',...
		];
end
%% Mean rate
mustr		=	['\t\t mu[i]\t<-   * ((1-' alphastrmu ') * (' fun '( ' par ' )) + ' alphastrmu ' ) \r\n'];

%% Response distribution
switch respDist
	case 'bernouilli'
		distr		= '\t\t y[i] ~ dbern(mu[i])\r\n';
	case 'binomial'
		distr		= '\t\t y[i,1] ~ dbin(mu[i],y[i,2])\r\n';
end

%% Hyperpriors
if ngroups>1
hyperprior = ['\t# threshold/lambda hyperprior \r\n',...
	'\t mulambda\t\t\t~ dnorm(0,1/(10)^2)\r\n',...
	'\t sigmalambda\t\t~ dunif(1.0E-3,1000)\r\n'];
else
hyperprior = ['\t# threshold/lambda hyperprior \r\n',...
	'\t mulambda\t\t\t<- 0\r\n',...
	'\t sigmalambda\t\t<- 1/(10)^2\r\n'];
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
	alphastr1,...
	'\t\t lambda[j]\t~ dnorm(mulambda,sigmalambda)\r\n',...
	'\t}\r\n',...
	'\t# constant/default width at F-1(ALPHA) and F-1(1-ALPHA) \r\n',...
	'\t alpha\t\t\t<- 0.1 \r\n',...
	alphastr2,...
	hyperalphastr,...
	hyperprior,...
	'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen(modelname,'w');
fprintf(fid,str);
fclose(fid);


keyboard

function p = psifun(x,lambda,alpha,varargin)
% P = PSIFUN(X,THETA,OMEGA,GAMMA,LAMBDA,ALPHA)
%
% Psychometric function of the form:
%
% p = (1-lambda)*( (1-alpha)*fun(x,lambda,omega,alpha)+alpha ) + alpha*lambda;
%
% See also PSIFIT, LOGISTICFUN, PROBITFUN, GUMBELFUN, REVGUMBELFUN,
% WEIBULLFUN, REVWEIBULLFUN

%% Initialization
hfun = keyval('function',varargin,@expfun);

if ~isa(hfun,'function_handle')
	error('Please pass a function handle')
end

	p = (1-lambda)*( (1-alpha)*hfun(x,lambda)+alpha ) + alpha*lambda;

