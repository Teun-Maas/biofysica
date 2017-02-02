% Jags-Ybinom-XnomSsubjCcat-MbinomBetaOmegaKappa.R
function samples = binorate_contrast(varargin)
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

%% Data
% 1: trial
% 2: condition, 1=visual left, 2=visual right, 3 = catch trial, 4 =
% visual left tactile right, 5  = visual right tactile right, 6 = tactile right
% 3: response 1=seen, 2 = unseen
% 4: correctness 1=correct, 0 = incorrect
% 5: reaction time (ms), determined via Matlab's clock
% 6: contrast 0 to 1, (Gabor-patch, Gaussian blur)
% 7: stimulus duration (ms)
% 8: 0
% 9: 0

cd('/Users/marcw/DATA/Mitchel Stokkermans');

d = dir('*.*.mat');
nfiles = numel(d);
z = [];
N = [];
C = [];
s = [];
for ii = 1:nfiles
	fname	= d(ii).name;
	load(fname,'Datamatrix'); % Datamatrix is interesting
	z	= [z; Datamatrix(3,:)']; %#ok<*NODEF>
	N	= [N; ones(size(Datamatrix(3,:)))'];
	s	= [s; repmat(ii,size(Datamatrix(3,:)))'];
	C	= [C; Datamatrix(2,:)'];
end
z = z==1;
z = double(z);

%% catch trials
% Determine 'guess' rate
sel			= ismember(C,[3 6]);
g			= z(sel); % Bernoulli = 0 or 1
Ng			= N(sel); % Bernoulli = 1 trial
sg			= s(sel);
[~,~,subs] = unique(sg);
g			= accumarray(subs,g,[],@sum); % Binomial: number of successes
Ng			= accumarray(subs,Ng,[],@sum); % Binomial: number of trials

%% Real trials
% [1 2]
% [4]
% [5]
% first redefine conditions
selVonly			= ismember(C,[1 2]);
selVT			= ismember(C,[4 5]);
selVlTr			= C==4;
selVrTr			= C==5;
selother		= ~ismember(C,[1 2 4 5]);
C(selVonly)		= 1;
C(selVlTr)		= 2;
C(selVrTr)		= 3;
C(selother)		= [];
z(selother)		= [];
N(selother)		= [];
s(selother)		= [];
selVT(selother) = [];
whos C z N s
sum(selVT)
C = [C; repmat(4,sum(selVT),1)];
z = [z; z(selVT)];
N = [N; N(selVT)];
s = [s; s(selVT)];


uS				= unique(s);
Nsubj			= numel(uS);

uC				= unique(C);
Ncond			= numel(uC);
zMatrix = NaN(Nsubj,Ncond);
nMatrix = NaN(Nsubj,Ncond);
sMatrix = NaN(Nsubj,Ncond);

for ii = 1:Ncond
	sel = C==uC(ii);
	ztmp			= z(sel);
	Ntmp			= N(sel);
	stmp			= s(sel);
	[stmp,~,subs]	= unique(stmp);
	zMatrix(:,ii)			= accumarray(subs,ztmp,[],@sum);
	nMatrix(:,ii)			= accumarray(subs,Ntmp,[],@sum);
	sMatrix(:,ii)			= stmp;
	cMatrix(:,ii)			= repmat(uC(ii),Nsubj,1);
	
end


% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,20000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

% graphic flags
diagFlag		= keyval('showDiag',varargin,true); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,true); % show posterior estimate distribution
predFlag		= keyval('showPred',varargin,true); % show posterior predictive
centroidFlag	= keyval('showCentroid',varargin,'mean'); % mode, median, mean


%% Actual regression
samples				= genMCMC(zMatrix,nMatrix,cMatrix,g,Ng,...
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

%%
close all
[Nsamples,Nsubj,Ncond] = size(samples.theta)

mu						= NaN(Nsubj,Ncond);
hdi						= NaN(Nsubj,Ncond,2);
mug						= NaN(Nsubj,1);
hdig						= NaN(Nsubj,2);

for ii = 1:Nsubj
	b		= samples.gamma(:,ii);
	mug(ii) = mean(b);
	hdig(ii,:)	= HDIofMCMC(b);
	
	for jj = 1:Ncond
		a			= squeeze(samples.theta(:,ii,jj));
		% 		a			= (1-b).*a+b;
		hdi(ii,jj,:)	= HDIofMCMC(a);
		mu(ii,jj)		= mean(a);
	end
end

%%
figure(100)
subplot(121)
[~,idx] = sort(mu(:,1));

x = 1:Nsubj;
x = x';
col = [1 1 1; 0.7 0.7 0.7;0.7 0.7 0.7; 0 0 0];
lw = [2 1 1 2];
lc = [0 0 0; 0.7 0.7 0.7;0.7 0.7 0.7;0 0 0];
for ii = 1:Ncond
	y			= mu(idx,ii);
	e(:,1)		= y-squeeze(hdi(idx,ii,1));
	e(:,2)		= squeeze(hdi(idx,ii,2))-y;
	figure(100)
	if ismember(ii,[1 4])
		he(ii) = errorbar(x+0.1*ii-Ncond*0.1/2-0.1,y,e(:,1),e(:,2),'ko-','MarkerFaceColor',col(ii,:),'LineWidth',lw(ii),'Color',lc(ii,:));
	else
		he(ii) = plot(x+0.1*ii-Ncond*0.1/2-0.1,y,'k-','LineWidth',lw(ii),'Color',lc(ii,:));
	end
	hold on
end
y			= mug(idx);
e(:,1)		= y-squeeze(hdig(idx,1));
e(:,2)		= squeeze(hdig(idx,2))-y;
figure(100)
he(5) = errorbar(x,y,e(:,1),e(:,2),'kd-','MarkerFaceColor','w');

axis square
box off
xlim([0.5 8.5]);
ylim([-0.05 1.05]);
set(gca,'TickDir','out','XTick',x,'XTickLabel',x,'YTick',0:0.1:1);
h = horline(0,'k-');
set(h,'Color',[.7 .7 .7])
h = horline(1,'k-');
set(h,'Color',[.7 .7 .7])
hold on
xlabel('Subject');
ylabel('\theta | P(Seen)');
str = {'V','V_LT_R','V_RT_R','VT','False Alarm | T'};
legend(he,str,'Location','E')


figure(200)
col = parula(Nsubj);
for ii = 1:Nsubj
	
	x			= squeeze(samples.theta(:,ii,1));
	y			= squeeze(samples.theta(:,ii,2));
	[m,s,a] = ellipse(x,y);
	subplot(131)
	plot(x(1:100:end),y(1:100:end),'.','Color',col(ii,:));
	hold on
	ellipseplot(m,1.96*s,a,'Color',col(ii,:));
	
	x			= squeeze(samples.theta(:,ii,1));
	y			= squeeze(samples.theta(:,ii,3));
	[m,s,a] = ellipse(x,y);
	subplot(132)
	plot(x(1:100:end),y(1:100:end),'.','Color',col(ii,:));
	hold on
	ellipseplot(m,1.96*s,a,'Color',col(ii,:));
	
	x			= squeeze(samples.theta(:,ii,2));
	y			= squeeze(samples.theta(:,ii,3));
	[m,s,a] = ellipse(x,y);
	subplot(133)
	plot(x(1:100:end),y(1:100:end),'.','Color',col(ii,:));
	hold on
	ellipseplot(m,1.96*s,a,'Color',col(ii,:));
end

for ii = 1:3
	subplot(1,3,ii)
	axis([-0.1 1.1 -0.1 1.1]);
	axis square;
	h = unityline('k-');
	set(h,'Color',[.7 .7 .7]);
	if ii == 1
		xlabel('\theta_V');
		ylabel('\theta_{V_LT_R}');
	elseif ii ==2
		xlabel('\theta_V');
		ylabel('\theta_{V_RT_R}');
	elseif ii ==3
		xlabel('\theta_{V_LT_R}');
		ylabel('\theta_{V_RT_R}');
	end
	box off
	set(gca,'TickDir','out','XTick',0:0.1:1,'YTick',0:0.1:1);
	h = horline(0,'k-');
	set(h,'Color',[.7 .7 .7])
	h = horline(1,'k-');
	set(h,'Color',[.7 .7 .7])
	h = verline(0,'k-');
	set(h,'Color',[.7 .7 .7])
	h = verline(1,'k-');
	set(h,'Color',[.7 .7 .7])
end


figure(100)
subplot(122)
x			= squeeze(samples.theta(:,:,1));
y			= squeeze(samples.theta(:,:,4));
X			= y-x;
rope = 	mean(mean(samples.gamma));

plotPost(X(:),'ROPE',[-rope rope]);
hold on
xlim([-0.6 0.6])
horline(0,'k-');
% 	verline(0,'k-');
xlabel('\theta_{VT} - \theta_{V}');
ylabel('P(\theta_{VT} - \theta_{V})');
set(gca,'TickDir','out','XTick',-0.5:0.1:0.5);
%%
keyboard

function [samples,stats] = genMCMC(z,N,c,g,Ng,varargin)
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
parameters		= {'theta','omega','kappa',...
	'gamma','omegag','kappag'};
writemodel;


[Nsubj,Ncond]	= size(z)

dataStruct = struct('z',z,'N',N,'Nsubj',Nsubj,'c',c,'Ncond',Ncond,'g',g,'Ng',Ng);

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
for ii = 1:nChains
	initsStruct(ii).theta			= repmat(0.5,Nsubj,Ncond);
	initsStruct(ii).gamma			= repmat(0.5,Nsubj,1);
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
	'model{ \r\n',...
	'\tfor ( s in 1:Nsubj ) { \r\n',...
	'\tfor ( c in 1:Ncond ) { \r\n',...
	'\t\t\tz[s,c] ~ dbin( (1-gamma[s])*theta[s,c] + gamma[s], N[s,c] ) \r\n',...
	'\t\t\ttheta[s,c] ~ dbeta( omega[c]*(kappa[c]-2)+1 , (1-omega[c])*(kappa[c]-2)+1 ) \r\n',...
	'\t\t} \r\n',...
	'\t\tg[s] ~ dbin(gamma[s],Ng[s]) \r\n',...
	'\tgamma[s] ~ dbeta(omegag*(kappag-2)+1 , (1-omegag)*(kappag-2)+1) # broad uniform \r\n',...
	'\t} \r\n',...
	'\tfor ( c in 1:Ncond ) { \r\n',...
	'\tomega[c] ~ dbeta(1,1) # broad uniform \r\n',...
	'\tkappa[c] <- kappaMinusTwo[c] + 2\r\n',...
	'\tkappaMinusTwo[c] ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10\r\n',...
	'\t\t} \r\n',...
	'\tomegag ~ dbeta(1,1) # broad uniform \r\n',...
	'\tkappag <- kappaMinusTwog + 2\r\n',...
	'\tkappaMinusTwog ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10\r\n',...
	'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);


