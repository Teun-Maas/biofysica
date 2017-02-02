function visualspeechrecognition_lr2016(varargin)
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
datadir			= '/Volumes/mbaudit4/Marc van Wanrooij/data/words';
cd(datadir);
load('spin');
% spin T C SNR L S

[uwords,~,words]	= unique(T);
% [uresponse,~,response]	= unique(R);

[~,~,subjects]		= unique(S);

% sel			= M==1;
sel = M==2;
y			= C(sel);
r			= R(sel);
t			= T(sel);
o			= O(sel);
s			= subjects(sel);
us			= unique(s);
ns			= numel(us);
w			= words(sel);
uw			= unique(w);
nw			= numel(uw);

z = NaN(ns,nw);
N = z;
g = z;
Ng = z;
for ii = 1:ns
	for jj = 1:nw
		selsub		= s==us(ii);
		selword		= w==uw(jj);
		sel			= selsub & selword;
		co			= unique(o(selword)); % current order
		selorder	= o==co;
		
		selresp		= strcmp(uwords(uw(jj)),r)';
		selstim		= strcmp(uwords(uw(jj)),t)';
		

		selfalse	= selsub & selresp & ~selstim & selorder';
	
		
		z(ii,jj)	= sum(y(sel));
		N(ii,jj)	= sum(sel);
		
		g(ii,jj) = sum(selfalse);
		Ng(ii,jj) = sum(selsub & ~selstim & selorder');
		% 		selfalse = s==us(ii) & strcmp(uwords(uw(jj)),r);
		% 		sum(selfalse)
		
% 		a = sum(selresp);
% 		
% 		b = sum(selstim);
% 		c = sum(selword);
% 		d = sum(~selstim);
% 		% e = sum(selword&selstim);
% 		e = sum(selresp & ~selstim);
% 		f = sum(~selstim & selorder');
% 		[a b c d e f]
	end
end

[~,idx]	= sort(mean(z./N,2));
zsorted		= z(idx,:);
Nsorted		= N(idx,:);
gsorted		= g(idx,:);
Ngsorted		= Ng(idx,:);

[~,idx]	= sort(mean(zsorted./Nsorted,1));
zsorted		= zsorted(:,idx);
Nsorted		= Nsorted(:,idx);
gsorted		= gsorted(:,idx);
Ngsorted		= Ngsorted(:,idx);

figure(1)
subplot(131)
imagesc(zsorted./Nsorted)
axis square
set(gca,'YDir','normal');
colorbar
caxis([0 1]);

subplot(132)
imagesc(gsorted./Ngsorted)
axis square
set(gca,'YDir','normal');
colorbar
caxis([0.05 0.15]);
% z = zsorted;
% N = Nsorted;
% g = gsorted;
% Ng = Ngsorted;

sum(gsorted./Ngsorted,2)



theta = z./N;
gamma = g./Ng;
theta = (theta-gamma)./(1-gamma);

subplot(133)
imagesc(theta)
axis square
set(gca,'YDir','normal');
colorbar
caxis([0 1]);


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


%% Actual regression
samples				= genMCMC(z,N,g,Ng,...
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
	pause

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
mutheta		= squeeze(mean(samples.theta));
mugamma		= squeeze(mean(samples.gamma));
mup			= squeeze(mean(samples.p));
muw			= squeeze(mean(samples.w));
rho			= z./N;

[~,idx]	= sort(mup);
mutheta		= mutheta(idx,:);
mugamma		= mugamma(idx,:);
rho			= rho(idx,:);

[~,idx]	= sort(muw);
mutheta		= mutheta(:,idx);
mugamma		= mugamma(:,idx);
rho			= rho(:,idx);

figure(666)
clf
subplot(131)
imagesc(rho)
hold on
axis square
set(gca,'YDir','normal','TickDir','out');
box off
caxis([0 1])
xlabel('Word');
ylabel('Subject');
title('Mean performance (%)');



subplot(132)
imagesc(mutheta)
hold on
contour(mutheta,0:0.2:1,'k')
xlabel('Word');
ylabel('Subject');
title('Posterior rate');

% plot(muw/max(muw)*18,'w')
% plot(mup/max(mup)*50,1:19,'w')

axis square
set(gca,'YDir','normal','TickDir','out');
box off;
caxis([0 1])

subplot(133)
imagesc(mugamma)
hold on
axis square
set(gca,'YDir','normal','TickDir','out');
box off
caxis([0 1])
xlabel('Word');
ylabel('Subject');
title('Posterior guess rate');

pa_datadir;
print('-dpng',[mfilename '_1']);

%% Separable rates for Person/subject and Words
figure(777)

p			= samples.p;
w			= samples.w;
np			= size(p,2);
nw			= size(w,2);

for ii = 1:np
	postSummaryP(ii) = summarizepost(p(:,ii));
end
mu = [postSummaryP.mean];
E = [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
[~,idx] = sort(mu);
mu		= mu(idx);
E = E(:,idx);
x = 1:np;

subplot(121)
cla
errorpatch(x,mu,E);
hold on
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'YTick',0:0.2:1,'TickDir','out','XTick',1:3:np);
box off
xlabel('Subject');
ylabel('Subjective visual recognition rate');
xlim([0 np+1]);
ylim([-0.1 1.1])
horline([0.1  1],'k:');
bf_text(0.05,0.95,'A');
plot([1 np],[0.1 1],'k:');

for ii = 1:nw
	postSummaryW(ii) = summarizepost(w(:,ii));
end

mu = [postSummaryW.mean];
E = [[postSummaryW.hdiLow]; [postSummaryW.hdiHigh]];
[~,idx] = sort(mu);
mu		= mu(idx);
E = E(:,idx);
subplot(122)
cla
errorpatch(1:nw,mu,E);
hold on
plot(1:nw,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'YTick',0:0.2:1,'TickDir','out')
words = unique(T);
words = words(idx);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);

box off
xlabel('Word');
ylabel('Visual word recognition rate');
xlim([0 nw+1]);
ylim([-0.1 1.1])
horline([0.1 1],'k:');
plot([1 50],[0.1 1],'k:');
bf_text(0.05,0.95,'B');

pa_datadir;
print('-dpng',[mfilename '_2']);

% Remaining questions:
% - if separable, does this mean, that we can suffice with measuring 1
% subject with all words, and many subjects with 1 word?
% - does order in sentence matter?
% - Why is it that some words are easier to recognize than others
% - is model correct?i.e. is gamma rate correctly used?
% - Subject characteristics? Audiograms? Age?
% - Use hierarchical priors?
% - Is data correct? i.e. does every subject have total set? Is there
% some data missing?

%%
keyboard


function [samples,stats] = genMCMC(z,N,g,Ng,varargin)
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
% fun				= keyval('function',varargin);
runjagsMethod	= keyval('runjagsMethod',varargin);
dic				= keyval('dic',varargin);
modelname		= fullfile(pwd, 'model.txt');

%% Write the model
% first check parameters to be monitored
parameters		= {'theta','p','w','gamma'};
writemodel;


% Specify data, as a structure
[np,nw] = size(z);
dataStruct = struct('z',z,'N',N,'np',np,'nw',nw,'g',g,'Ng',Ng);

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
	initsStruct(ii).gamma				= repmat(0.1,np,nw);  % ~std
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
% 	save([saveName 'Mcmc'],'samples');
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
	'\t    g[i,j] ~ dbin(gamma[i,j],Ng[i,j]) \r\n',...
	'\t# Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t    theta[i,j] <- (1-gamma[i,j])*p[i]*w[j] +gamma[i,j]\r\n',...
	'\t    gamma[i,j] ~ dbeta(1,1) \r\n',...
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


