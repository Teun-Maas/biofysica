function samples = dexppsifit(x,y,s,e,varargin)
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
	e = e';
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

	samples = jags_expfit(x,y,s,e,...
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




function [samples,stats] = jags_expfit(x,y,s,e,varargin)
% SAMPLES = JAGS_PSIFIT(X,Y)
%
% See also PSIFIT

%% Psychometric function
fun				= keyval('function',varargin,@expfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
alpha			= keyval('alpha',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 alpha for all subjects/conditions
% modelname		= fullfile(pwd, 'jags_exp_model.txt');
modelname = '/Users/marcw/Dropbox/Manuscript/Luuk van de Rijt/#4 SPIN AV CI/matlab/dexp_model.txt';
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
parameters		= {'lambda','dlambda','mulambda','mudlambda','alpha','mualpha'};
dataStruct		= struct('x',x,'y',y,'s',s,'e',e,'Ntotal',Ntotal,'Nsubj',Ngroups); % data structure
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




