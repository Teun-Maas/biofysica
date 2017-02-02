function auditoryspeechrecognition_lr2016(varargin)
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


loadFlag		= keyval('load',varargin,false);
sampleFlag		= keyval('sample',varargin,true);
datadir			= '/Volumes/mbaudit4/Marc van Wanrooij/data/words';
cd(datadir);

if ~loadFlag
	load('spin');
	[uwords,~,words]	= unique(T);
	[~,~,subjects]		= unique(S);
	
	sel			= M==1;
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
	snr			= SNR(sel);
	usnr		= unique(snr);
	nsnr		= numel(usnr);
	
	%% Initialization
	z			= NaN(ns,nw,nsnr);
	N			= z;
	g			= z;
	Ng			= z;
	xSNR		= z;
	xS			= z;
	xW			= z;
	
	%% Loop
	for ii = 1:ns % for every subject
		for jj = 1:nw % for every word
			for kk = 1:nsnr % for every signal-to-noise ratio
				selsub		= s==us(ii);
				selword		= w==uw(jj);
				selsnr		= snr==usnr(kk);
				sel			= selsub & selword & selsnr';
				
				co			= unique(o(selword)); % current order
				selorder	= o==co;
				selresp		= strcmp(uwords(uw(jj)),r)';
				selstim		= strcmp(uwords(uw(jj)),t)';
				selfalse	= selsub & selresp & ~selstim & selorder' &selsnr';
				
				
				z(ii,jj,kk)		= sum(y(sel));
				N(ii,jj,kk)		= sum(sel);
				g(ii,jj,kk)		= sum(selfalse);
				Ng(ii,jj,kk)	= sum(selsub & ~selstim & selorder' & selsnr');
				xSNR(ii,jj,kk)	= usnr(kk);
				xS(ii,jj,kk)	= ii;
				xW(ii,jj,kk)	= jj;
			end
		end
	end
	%% Plot raw data
	close all
	figure(1)
	rho = z./N;
	for ii = 1:nsnr
		r = squeeze(rho(:,:,ii));
		squeeze(N(:,:,ii))
		
		subplot(2,3,ii)
		imagesc(r);
		box off
		axis square
		caxis([0 1]);
		set(gca,'YDir','normal','TickDir','out');
		
		subplot(2,3,6)
		m = mean(r);
		plot(usnr(ii),m(:),'.');
		hold on
		
		m = mean(r,2);
		plot(usnr(ii),m(:),'.');
		hold on
		xlim([-25 0]);
	end
	
	z		= z(:);
	N		= N(:);
	g		= g(:);
	Ng		= Ng(:);
	xSNR	= xSNR(:);
	xS		= xS(:);
	xW		= xW(:);
	
	sel		= N~=0;
	z		= z(sel);
	N		= N(sel);
	g		= g(sel);
	Ng		= Ng(sel);
	xSNR	= xSNR(sel);
	xS		= xS(sel);
	xW		= xW(sel);
	
	data.hit			= z;
	data.ntrials		= N;
	data.catchhit		= g;
	data.ncatchtrials	= Ng;
	data.subject		= xS;
	data.word			= xW;
	data.SNR			= xSNR; %#ok<STRNU>
	
	save('spinauditoryraw','data');
elseif loadFlag
		load('spin');

	load('spinauditoryraw');
	z = data.hit; %#ok<NODEF>
	N = data.ntrials;
	g = data.catchhit;
	Ng = data.ncatchtrials;
	xS = data.subject;
	xW = data.word;
	xSNR = data.SNR;
end

if sampleFlag
	%% MCMC
	chain_globals;
	numSavedSteps	= keyval('numSavedSteps',varargin,1000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
	thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
	burnInSteps		= keyval('burnInSteps',varargin,100);
	saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
	nChains			= keyval('nChains',varargin,nChainsDefault);
	runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
	dic				= keyval('dic',varargin,false);
	
	% graphic flags
	diagFlag		= keyval('showDiag',varargin,true); % show MCMC diagnostics
	postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
	
	
	%% Actual regression
	samples				= genMCMC(z,N,g,Ng,xS,xW,xSNR,...
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
		plotdiag(samples);
	end
	
	%% Extract chain values:
	samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D
	
	save('spinauditoryMCMC','samples');
	
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
	
elseif ~sampleFlag
	
	load('spinauditoryMCMC');
end

%%

% keyboard
%% Separable rates for Person/subject and Words
figure(777)

% Theta
p			= samples.ptheta;
w			= samples.wtheta;
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

subplot(231)
cla
errorpatch(x,mu,E);
hold on
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out','XTick',1:3:np,'XTickLabel',idx(1:3:np));
box off
xlabel('Subject');
ylabel('Subjective auditory threshold (dB)');
xlim([0 np+1]);
ylim([-15 5])
% horline([0.1  1],'k:');
bf_text(0.05,0.95,'A');
% plot([1 np],[0.1 1],'k:');

for ii = 1:nw
	postSummaryW(ii) = summarizepost(w(:,ii));
end

mu = [postSummaryW.mean];
E = [[postSummaryW.hdiLow]; [postSummaryW.hdiHigh]];
[~,idx] = sort(mu);
mu		= mu(idx);
E = E(:,idx);
subplot(232)
cla
errorpatch(1:nw,mu,E);
hold on
plot(1:nw,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out')
words = unique(T);
words = words(idx);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);

box off
xlabel('Word');
ylabel('Word auditory threshold (dB)');
xlim([0 nw+1]);
ylim([-15 5])
% horline([0.1 1],'k:');
% plot([1 50],[0.1 1],'k:');
bf_text(0.05,0.95,'B');


% Omega
p			= samples.pomega;
w			= samples.womega;
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

subplot(234)
cla
errorpatch(x,mu,E);
hold on
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out','XTick',1:3:np,'XTickLabel',idx(1:3:np));
box off
xlabel('Subject');
ylabel('Subjective auditory width (dB)');
xlim([0 np+1]);
ylim([-2 10])
% horline([0.1  1],'k:');
bf_text(0.05,0.95,'C');
% plot([1 np],[0.1 1],'k:');

for ii = 1:nw
	postSummaryW(ii) = summarizepost(w(:,ii));
end

mu = [postSummaryW.mean];
E = [[postSummaryW.hdiLow]; [postSummaryW.hdiHigh]];
[~,idx] = sort(mu);
mu		= mu(idx);
E = E(:,idx);
subplot(235)
cla
errorpatch(1:nw,mu,E);
hold on
plot(1:nw,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out')
words = unique(T);
words = words(idx);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);

box off
xlabel('Word');
ylabel('Word auditory width (dB)');
xlim([0 nw+1]);
ylim([-2 10])
% horline([0.1 1],'k:');
% plot([1 50],[0.1 1],'k:');
bf_text(0.05,0.95,'D');


% Gamma
g			= samples.gamma;

[~,np,nw]			= size(g);

mu	= NaN(np,nw);
L	= NaN(np,nw);
U	= NaN(np,nw);
for ii = 1:np
	for jj = 1:nw
		postSummaryG(ii,jj) = summarizepost(squeeze(g(:,ii,jj)));
		mu(ii,jj) = postSummaryG(ii,jj).mean;
		L(ii,jj) = postSummaryG(ii,jj).hdiLow;
		U(ii,jj) = postSummaryG(ii,jj).hdiHigh;
		
	end
end


[~,idx1] = sort(mean(mu));
mu		= mu(:,idx1);
[~,idx2] = sort(mean(mu,2));
mu		= mu(idx2,:);
subplot(236)
imagesc(mu);
axis square;
colorbar;
caxis([0 0.2])
set(gca,'TickDir','out','YDir','normal');
box off

load('spin');

words = unique(T);
words = words(idx1);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);
set(gca,'YTick',1:3:np,'YTickLabel',idx2(1:3:np));



% pa_datadir;
% print('-dpng',[mfilename '_2']);

%%
% To do:
% Posterior predictive

%%
keyboard


function [samples,stats] = genMCMC(z,N,g,Ng,xS,xW,xSNR,varargin)
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
parameters		= {'ptheta','wtheta','pomega','womega'};
writemodel;


% Specify data, as a structure
np = max(xS);
nw = max(xW);

dataStruct = struct('z',z,'N',N,'np',np,'nw',nw,'x',xSNR,'p',xS,'w',xW,'nTotal',numel(z));

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
for ii = 1:nChains
	initsStruct(ii).ptheta				= zeros(np,1);
	initsStruct(ii).wtheta				= zeros(nw,1);
	initsStruct(ii).pomega				= repmat(0.5,np,1);  % ~std
	initsStruct(ii).womega				= repmat(0.5,nw,1);  % ~std
	initsStruct(ii).gamma				= repmat(0.1,np,nw);  % ~std
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


parameterNames		= fieldnames(samples) % get all parameter names


%% Save samples
if ~isempty(saveName)
	% 	save([saveName 'Mcmc'],'samples');
end


function writemodel
%%

str = [
	'\tmodel{ \r\n',...
	'\tfor (i in 1:nTotal){ \r\n',...
	'\t# Word performance Is Binomially Distributed \r\n',...
	'\t z[i] ~ dbin(theta[i],N[i]) \r\n',...
	'\t# Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t theta[i]  <- (1-0.1)*f[i] +0.1\r\n',...
	'\t	f[i]	<- ilogit( (2*log(1/alpha-1))/(pomega[p[i]]+womega[w[i]]) * (x[i]-ptheta[p[i]]-wtheta[w[i]]) )  \r\n',...
	'\t} \r\n',...
	'\t# Priors For People and Words \r\n',...
	'\tfor (i in 1:np){ \r\n',...
	'\t\t pomega[i]\t~ dgamma(1,1)  \r\n',...
	'\t\t ptheta[i]\t~ dnorm(0,1/(10)^2)\r\n',...
	'\t} \r\n',...
	'\tfor (j in 1:nw){ \r\n',...
	'\t\t womega[j]\t~ dgamma(1,1)  \r\n',...
	'\t\t wtheta[j]\t~ dnorm(0,1/(10)^2)\r\n',...
	'\t} \r\n',...
	'\t# constant/default width at F-1(ALPHA) and F-1(1-ALPHA) \r\n',...
	'\t alpha\t\t\t<- 0.1 \r\n',...
	'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);

function plotdiag(samples)
parameterNames	= fieldnames(samples); % get all parameter names
for parIdx			= 1:numel(parameterNames)
	n = size(samples.(parameterNames{parIdx}));
	if(length(n))<4
		n = size(samples.(parameterNames{parIdx}),3);
		
		for ii = 1:2
			figure
			a		= squeeze(samples.(parameterNames{parIdx})(:,:,ii));
			samp	= samples;
			samp.(parameterNames{parIdx}) = a;
			
			diagmcmc(samp,'parName',parameterNames{parIdx});
		end
	end
end
