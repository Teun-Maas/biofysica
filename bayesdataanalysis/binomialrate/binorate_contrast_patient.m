function samples = binorate_contrast_patient(varargin)
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

cd('/Users/marcw/DATA/Ahmed Gardoh');
cd('/Users/marcw/DATA/Mitchel Stokkermans');
% cd('/Users/marcw/DATA/Ahmed Gardoh/Ahmed copy');
cd('/Users/marcw/DATA/Ahmed Gardoh/Healthy people');

d = dir('*.mat');
nfiles = numel(d);
z = [];
N = [];
C = [];
s = [];
dT = [];
for ii = 1:nfiles
	fname	= d(ii).name;
	load(fname,'Datamatrix'); % Datamatrix is interesting
	sel = ismember(Datamatrix(2,:),1:6);
	Datamatrix = Datamatrix(:,sel);
	z	= [z; Datamatrix(3,:)']; %#ok<*NODEF>
	N	= [N; ones(size(Datamatrix(3,:)))'];
	s	= [s; repmat(ii,size(Datamatrix(3,:)))'];
	C	= [C; Datamatrix(2,:)'];
	dT = [dT; Datamatrix(10,:)'];
end
z = z==1;
z = double(z);

% unique(dT)
%%
sel = dT>0;
sum(sel)

% z = z(sel);
% N = N(sel);
% C = C(sel);
% s = s(sel);

%% catch trials
% Determine 'guess' rate
sel			= ismember(C,3);
g			= z(sel); % Bernoulli = 0 or 1
Ng			= N(sel); % Bernoulli = 1 trial
sg			= s(sel);
[~,~,subs] = unique(sg);
g			= accumarray(subs,g,[],@sum); % Binomial: number of successes
Ng			= accumarray(subs,Ng,[],@sum); % Binomial: number of trials

%% Real trials
uS				= unique(s);
Nsubj			= numel(uS);

uC				= unique(C);
Ncond			= numel(uC)
zMatrix			= NaN(Nsubj,Ncond);
nMatrix			= NaN(Nsubj,Ncond);
sMatrix			= NaN(Nsubj,Ncond);

for ii = 1:Ncond
	sel					= C==uC(ii);
	ztmp				= z(sel);
	Ntmp				= N(sel);
	stmp				= s(sel);
	[stmp,~,subs]		= unique(stmp);
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
% diagFlag		= keyval('showDiag',varargin,true); % show MCMC diagnostics
% postFlag		= keyval('showPost',varargin,true); % show posterior estimate distribution
% predFlag		= keyval('showPred',varargin,true); % show posterior predictive
% centroidFlag	= keyval('showCentroid',varargin,'mean'); % mode, median, mean


%% Actual regression
samples				= genMCMC(zMatrix,nMatrix,cMatrix,g,Ng,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic);



%% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

%%
close all
[~,Nsubj,Ncond] = size(samples.theta);

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
		hdi(ii,jj,:)	= hdimcmc(a);
		mu(ii,jj)		= mean(a);
	end
end


%%
[~,~,idx] = unique(mu(:,1));
idx = zeros(size(idx));
% close all
% % plot(mu(:,1),idx)
% plot(mu(idx,1))

% idx = 40*mu(:,1);
close all
figure(101)
cnt = 0;
clear x y e h
cndIdx = [1 4];
subplot(152)
for sIdx = 1:Nsubj
	cnt = cnt+1;
	x			= 1:2;
	y			= mu(sIdx,cndIdx);
	e(:,1)		= y-squeeze(hdi(sIdx,cndIdx,1));
	e(:,2)		= squeeze(hdi(sIdx,cndIdx,2))-y;
	h(cnt) = errorbar(x+0.01*idx(sIdx),y,e(:,1),e(:,2),'o-','MarkerFaceColor','w'); %#ok<*AGROW>
	hold on
end
ylim([-0.1 1.1])
xlim([0.5 2.5])
set(gca,'TickDir','out',...
	'XTick',1:2,'XTickLabel',{'V_L','V_LT_L'},...
	'YTick',0:0.2:1,'YTickLabel',[]);
box off
horline([0 1]);
xlabel('Condition');
% ylabel('Visual recognition rate \theta');
bf_text(0.1,0.95,char(66));
title('Congruent');

cnt = 0;
clear x y e h
cndIdx = [2 5];
subplot(153)
for sIdx = 1:Nsubj
	cnt = cnt+1;
	x			= 1:2;
	y			= mu(sIdx,cndIdx);
	e(:,1)		= y-squeeze(hdi(sIdx,cndIdx,1));
	e(:,2)		= squeeze(hdi(sIdx,cndIdx,2))-y;
	h(cnt) = errorbar(x+0.01*idx(sIdx),y,e(:,1),e(:,2),'o-','MarkerFaceColor','w');
	hold on
end
ylim([-0.1 1.1])
xlim([0.5 2.5])
set(gca,'TickDir','out',...
	'XTick',1:2,'XTickLabel',{'V_R','V_RT_L'},...
	'YTick',0:0.2:1,'YTickLabel',[]);
box off
horline([0 1]);
xlabel('Condition');
% ylabel('Visual recognition rate \theta');
bf_text(0.1,0.95,char(67));
title('Incongruent');


cnt = 0;
clear x y e h
cndIdx = [3 6];
subplot(151)
for sIdx = 1:Nsubj
	cnt = cnt+1;
	x			= 1:2;
	y			= mu(sIdx,cndIdx);
	e(:,1)		= y-squeeze(hdi(sIdx,cndIdx,1));
	e(:,2)		= squeeze(hdi(sIdx,cndIdx,2))-y;
	h(cnt) = errorbar(x+0.01*idx(sIdx),y,e(:,1),e(:,2),'o-','MarkerFaceColor','w');
	hold on
end
ylim([-0.1 1.1])
xlim([0.5 2.5])
set(gca,'TickDir','out',...
	'XTick',1:2,'XTickLabel',{'False alarm','T_L'},...
	'YTick',0:0.2:1);
box off
horline([0 1]);
xlabel('Condition');
ylabel('Visual recognition rate \theta');
legend(h,num2str((1:Nsubj)'),'Location','NE');
bf_text(0.1,0.95,char(65));
title({'False alarm','Distractor'});

cnt = 0;
clear x y e h
cndIdx = [1 2];
subplot(154)
for sIdx = 1:Nsubj
	cnt = cnt+1;
	x			= 1:2;
	y			= mu(sIdx,cndIdx);
	e(:,1)		= y-squeeze(hdi(sIdx,cndIdx,1));
	e(:,2)		= squeeze(hdi(sIdx,cndIdx,2))-y;
	h(cnt) = errorbar(x+0.01*idx(sIdx),y,e(:,1),e(:,2),'o-','MarkerFaceColor','w');
	hold on
end
ylim([-0.1 1.1])
xlim([0.5 2.5])
set(gca,'TickDir','out',...
	'XTick',1:2,'XTickLabel',{'V_L','V_R'},...
	'YTick',0:0.2:1,'YTickLabel',[]);
box off
horline([0 1]);
xlabel('Condition');
% ylabel('Visual recognition rate \theta');
bf_text(0.1,0.95,char(68));
title('Side Visual');

cnt = 0;
clear x y e h
cndIdx = [4 5];
subplot(155)
for sIdx = 1:Nsubj
	cnt = cnt+1;
	x			= 1:2;
	y			= mu(sIdx,cndIdx);
	e(:,1)		= y-squeeze(hdi(sIdx,cndIdx,1));
	e(:,2)		= squeeze(hdi(sIdx,cndIdx,2))-y;
	h(cnt) = errorbar(x+0.01*idx(sIdx),y,e(:,1),e(:,2),'o-','MarkerFaceColor','w');
	hold on
end
ylim([-0.1 1.1])
xlim([0.5 2.5])
set(gca,'TickDir','out',...
	'XTick',1:2,'XTickLabel',{'V_LT_L','V_RT_L'},...
	'YTick',0:0.2:1,'YTickLabel',[]);
box off
horline([0 1]);
xlabel('Condition');
% ylabel('Visual recognition rate \theta');
bf_text(0.1,0.95,char(69));
title('Side Visuo-Tactile');

print('-depsc',[mfilename '_comparison']);

%%
figure(100)
clf
m	= NaN(Nsubj*2,2);
sd	= m;
a	= NaN(Nsubj*2,1);
hd	= m;
col = lines(Nsubj*2);
X = [];
Y = [];
for ii = [0 1]
	for sIdx = 1:Nsubj
		x = squeeze(samples.theta(:,sIdx,1+ii));
		y = squeeze(samples.theta(:,sIdx,4+ii));
		X = [X x];
				Y = [Y y];

		[m(sIdx+Nsubj*ii,:),sd(sIdx+Nsubj*ii,:),a(sIdx+Nsubj*ii)] = ellipse(x,y-x);
		hd(sIdx+Nsubj*ii,:) = hdimcmc(y-x);
	end
end
[~,idx] = sort(m(:,1)); % sort by strength
for ii = [0 1]
	for sIdx = 1:Nsubj
		x = X(:,idx(sIdx+Nsubj*ii));
		y = Y(:,idx(sIdx+Nsubj*ii));
		
% 		plot(x(1:10:end),y(1:10:end)-x(1:10:end),'.','Color',col(idx(sIdx+Nsubj*ii),:));
		hold on
		ellipseplot(m(idx(sIdx+Nsubj*ii),:),sd(idx(sIdx+Nsubj*ii),:),a(idx(sIdx+Nsubj*ii)),'Color',col(idx(sIdx+Nsubj*ii),:));
	end
end
xi = linspace(0,1,100);
N = 2;
P = polyfit(m(idx,1),m(idx,2),N);
yi = polyval(P,xi);
plot(xi,yi,'k-','LineWidth',2)
% errorbar(m(idx,1),m(idx,2),m(idx,2)-hd(idx,1),hd(idx,2)-m(idx,2),'ko','MarkerFaceColor','w');
plot(m(idx,1),m(idx,2),'ko','MarkerFaceColor','w');

axis square
xlim([-0.1 1.1]);
ylim([-0.1 0.4]);
horline;
xlabel('\theta_{V_L}')
ylabel('\Delta \theta (V_LT_L-V_L)');
set(gca,'TickDir','out',...
	'XTick',0:0.2:1,...
	'YTick',0:0.1:0.3);
box off;
title('Inverse Effectiveness');
plot([0 1],[1 0],'k--');
h = text(0.8,0.2,'Ceiling'); set(h,'Rotation',-63,'VerticalAlignment','Bottom');
h = text(0.5,-0.001,'No effect'); set(h,'VerticalAlignment','Top','HorizontalAlignment','center');
verline([0 1]);

print('-depsc',[mfilename '_inverse']);

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
nChains			= keyval('nChains',varargin);
runjagsMethod	= keyval('runjagsMethod',varargin);
dic				= keyval('dic',varargin);
modelname		= fullfile(pwd, 'model.txt');

%% Write the model
% first check parameters to be monitored
parameters		= {'theta','omega','kappa',...
	'gamma','omegag','kappag'};
writemodel;


[Nsubj,Ncond]	= size(z);

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


