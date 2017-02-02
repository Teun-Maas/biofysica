% Jags-Ybinom-XnomSsubjCcat-MbinomBetaOmegaKappa.R 
function samples = binorate(z,n,s,varargin)
% MCMC = BINORATE(Z,N,S)
%
% Determine binomial rate
%
% Z = number of successes
% N = number of trials
% S = subject ID (or condition)
%
% MCMC structure contains fields with MCMC samples:
% - theta: rate
%
% MCMC = BINORATE(K,N,S,'NAME',VALUE)
% Additional name-value pair inputs include:
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
%
% Accompanies the book:
%   Kruschke, J. K. (2014). Doing Bayesian Data Analysis: 
%   A Tutorial with R, JAGS, and Stan. 2nd Edition. Academic Press / Elsevier.
% source('DBDA2E-utilities.R')


%% Initialization
close all
datadir = '/Volumes/mbaudit1/Marc van Wanrooij/SPIN/words';
cd(datadir);
load('spin');
% spin T C SNR L S

[~,~,words] = unique(T);
[~,~,subjects] = unique(S);

sel			= M==2;
y			= C(sel);
s			= subjects(sel);
us			= unique(s);
ns			= numel(us);
w		= words(sel);
uw			= unique(w);
nw			= numel(uw);

z = NaN(ns,nw);
N = z;
for ii = 1:ns
	for jj = 1:nw
		sel = s==us(ii) & w==uw(jj);
		z(ii,jj) = sum(y(sel));
		N(ii,jj) = sum(sel);
	end
end

[mu,idx]	= sort(mean(z./N,2));
zsorted		= z(idx,:);
Nsorted		= N(idx,:);

[mu,idx]	= sort(mean(zsorted./N,1));
zsorted		= zsorted(:,idx);
Nsorted		= Nsorted(:,idx);


figure
imagesc(zsorted./Nsorted)
axis square
set(gca,'YDir','normal');

z = zsorted;
N = Nsorted;


%%


% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,10000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,3000);
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
samples				= genMCMC(z,N,...
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
	parameterNames	= fieldnames(samples); % get all parameter names

	n = numel(parameterNames);
	for parIdx			= 1:n
		a = samples.(parameterNames{parIdx});
		a = a(:);
		subplot(n,n,n*(parIdx-1)+parIdx);
		cla
		plotpost(a);
		
	end
end

%%
close all
mutheta = squeeze(mean(samples.theta));
mup = squeeze(mean(samples.p));
muw = squeeze(mean(samples.w));

figure(666)
clf
subplot(221)
imagesc(z./N)
hold on
contour(z./N,0:0.2:1,'k')
axis square
set(gca,'YDir','normal');
caxis([0 1])

subplot(222)
plot(mup,1:length(mup))
axis square
box off

subplot(223)
plot(muw)
axis square
box off

subplot(224)
imagesc(mutheta)
hold on
contour(mutheta,0:0.2:1,'k')
axis square
set(gca,'YDir','normal');
caxis([0 1])


figure
imagesc(mutheta-z./N)
hold on
contour(mutheta-z./N,-1:1:1,'k')
colormap('jet')
axis square
set(gca,'YDir','normal');
colorbar
caxis([-1 1])

function [samples,stats] = genMCMC(z,N,varargin)
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
modelname		= fullfile(pwd, 'model.txt');

%% Write the model
% first check parameters to be monitored
	parameters		= {'theta','p','w'};
writemodel;


% Specify data, as a structure
[np,nw] = size(z)
dataStruct = struct('z',z,'N',N,'np',np,'nw',nw);

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
for ii = 1:nChains
% 	initsStruct(ii).theta			= zeros(np,nw); 
	initsStruct(ii).p				= repmat(0.5,np,1);  % ~std
	initsStruct(ii).w				= repmat(0.5,nw,1);  % ~std
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
if any(strcmp(func2str(fun),{'weibullfun','revweibullfun'}));
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


for ii = 1:ns % for every subject/group
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
	
	figure(200+figcnt);
	subplot(sb,sb,sbcnt)
	hold on
	
	%% Posterior credible psychometric curves
	halfway = ((1-clambda)+cgamma)/2;
	for cIdx	= cVec
		ypred	= psifun(xComb,theta(cIdx,ii),omega(cIdx,ii),gamma(cIdx,ii),lambda(cIdx,ii),0.1,'function',fun);
		xInt	= theta(cIdx,ii);

		plot(xComb , ypred,'k-','LineWidth',1.5,'Color',[.7 .7 .7]);
		plot( [xInt xInt],[halfway(ii) -0.1], 'k:','Color',[.7 .7 .7]);
	end
	
	%%
	% 	keyboard
	
	%% Predictive posterior max?
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
	
	plot(ux,rprop,'ks','MarkerFaceColor','w','MarkerSize',5);
	
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

function writemodel
% Placeholder function to generate JAGS model

%% Model
str = [
'\tmodel{ \r\n',...
  '\t# Word performance Is Binomially Distributed \r\n',...
  '\tfor (i in 1:np){ \r\n',...
  '\t  for (j in 1:nw){ \r\n',...
  '\t    z[i,j] ~ dbin(theta[i,j],N[i,j]) \r\n',...
  '\t  } \r\n',...
  '\t} \r\n',...
  '\t# Probability Correct Is Product Of Word By Person Rates \r\n',...
  '\tfor (i in 1:np){ \r\n',...
  '\t  for (j in 1:nw){ \r\n',...
  '\t    theta[i,j] <- p[i]*w[j] \r\n',...
  '\t  } \r\n',...
  '\t} \r\n',...
  '\t# Priors For People and Words \r\n',...
  '\tfor (i in 1:np){ \r\n',...
  '\t  p[i] ~ dbeta(1,1) \r\n',...
  '\t} \r\n',...
  '\tfor (j in 1:nw){ \r\n',...
  '\t  w[j] ~ dbeta(1,1) \r\n',...
  '\t} \r\n',...
'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);


