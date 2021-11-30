function samples = psifit(x,y,s,varargin)
% MCMC = PSIFIT(X,Y,S)
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
% predFlag		= keyval('showPred',varargin,false); % show posterior predictive



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


fun				= keyval('function',varargin,@logisticfun); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun
fig				= keyval('figure',varargin,200); % logistic, probit, gumbell, inversegumbell,  weibull, inverseweibull + fun

gamma			= keyval('gamma',varargin,'infer'); % 0.5 = 2AFC, 'infer' = estimate from data, 'single' = infer 1 gamma for all subjects/conditions
lambda			= keyval('lambda',varargin,'infer');% 0.001 = low lapse rate, 'infer' = estimate from data, 'single' = infer 1 lambda for all subjects/conditions

% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);
if exist('stan.m','file')
	def = 'stan';
else
	def='jags';
end
mcmcMethod		= keyval('mcmcMethod',varargin,def);

% graphic flags
diagFlag		= keyval('showDiag',varargin,false); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
predFlag		= keyval('showPred',varargin,false); % show posterior predictive
centroidFlag	= keyval('showCentroid',varargin,'mean'); % mode, median, mean

%% Actual regression

if strcmp(mcmcMethod,'stan')
	samples = stan_psifit(x,y,s,...
		'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
		'saveName',saveName,'nChains',nChains,'function',fun,...
		'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
		'dic',dic,'gamma',gamma,'lambda',lambda,'respDist',respDist);
else
	samples = jags_psifit(x,y,s,...
		'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
		'saveName',saveName,'nChains',nChains,'function',fun,...
		'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
		'dic',dic,'gamma',gamma,'lambda',lambda,'respDist',respDist);
end


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
		[r,ux]	= assemble(y(sel),x(sel),'fun',@sum);
		n		= assemble(ones(size(y(sel))),x(sel),'fun',@sum);
	elseif size(y,2)==2 % y = [rate n]
		 % % this should work
% 		r = y(sel,1);
% 		n = y(sel,2);
% 		ux = x(sel);
% % but perhaps the experimenter did not 'assemble' the data correctly
xs = x(sel);
xs = round(xs/2.5)*2.5;
		[r,ux]	= assemble(y(sel,1),xs,'fun',@sum);
		n		= assemble(y(sel,2),xs,'fun',@sum);		
	end
	rprop			= r./n;
	[lb, ub]		= binomialci(r, n, 0.05);	
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

