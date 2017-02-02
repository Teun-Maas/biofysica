function auditoryvisualspeechrecognition_lr2016(varargin)
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
[uwords,~,words]	= unique(T);
[usubjects,~,subjects]		= unique(S);
mod = 1;
[zA,NA,gA,NgA,xSNRA,xSA,xWA] = getdata(M,C,R,T,O,SNR,subjects,words,uwords,mod);

%%
usubjects = unique(xSA);
uwords = unique(xWA);
nsubjects = numel(usubjects);
nwords = numel(uwords);
% nsubjects = 2;
% nwords = 1;
close all
for subjectIdx = 5
for wordIdx = 1
% 	sel = xSA==usubjects(subjectIdx) & xWA==uwords(wordIdx);
		sel = xWA==uwords(wordIdx);

	sum(NA(sel))
	samples = psifit(xSNRA(sel),[zA(sel) NA(sel)],xSA(sel),'showDiag','true','gamma',0,'lambda',0);
end
end
%%
keyboard

mod = 2;
[zV,NV,gV,NgV,~,xSV,xWV] = getdata(M,C,R,T,O,SNR,subjects,words,uwords,mod);

mod = 3;
[zAV,NAV,gAV,NgAV,xSNRAV,xSAV,xWAV] = getdata(M,C,R,T,O,SNR,subjects,words,uwords,mod);

%% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,2000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,100);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

%% graphic flags
diagFlag		= keyval('showDiag',varargin,false); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution

%% Actual regression
samples				= genMCMC(zA,NA,gA,NgA,xSA,xWA,xSNRA,...
	zV,NV,gV,NgV,xSV,xWV,...
	zAV,NAV,gAV,NgAV,xSAV,xWAV,xSNRAV,...
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


%% Posterior estimates
if postFlag
	plotthispost(samples);
end

%% Separable rates for Person/subject and Words
plotparams(samples,T)
%%
% To do:
%% Posterior predictive
x = -22:2;
for jj = 5
	for ii = 1:16
		sub		= ii;
		word	= jj;
		[sub word]
		% auditory
		nsamples = size(samples.pthetaA,1);
		sel		= xSA==sub & xWA==word;
		strA			= ['z_A = ' num2str(zA(sel)') ', N_A = ' num2str(NA(sel)') ,...
			', z_\gamma_A = ' num2str(sum(gA(sel))') ', N_\gamma_A = ' num2str(sum(NgA(sel))') ];

		yA		= zA(sel)./NA(sel);
		xA		= xSNRA(sel);
		nA		= NA(sel);
		nx		= numel(xA);
		z		= NaN(nsamples,nx);
		for kk = 1:nsamples
			pthetaA			= samples.pthetaA(kk,sub);
			wthetaA			= samples.wthetaA(kk,word);
			pomegaA			= samples.pomegaA(kk,sub);
			womegaA			= samples.womegaA(kk,word);
			gammaA			= squeeze(samples.gammaA(kk,sub,word));
			
			theta			= logisticfun(xA,pthetaA+wthetaA,pomegaA+womegaA,0.1);
			theta			= gammaA+(1-gammaA).*theta;
			tmp				= binornd(nA,theta)./nA;
			z(kk,:)			= tmp;
		end
		X		= xA;
		u		= unique(z);
		[X,Y]	= meshgrid(X,u);
		Z		= NaN(size(X));
		for kk	= 1:numel(xA)
			for ll = 1:numel(u)
				sel =  z(:,kk)==u(ll);
				Z(ll,kk)	= round(sum(sel)*1000);
			end
		end
		Z = Z./max(Z(:));
		Z = Z*500;
		X = X(:);
		Y = Y(:);
		Z = Z(:);
		sel = isnan(Z) | Z<=0.9;
		sel = ~sel;
		X = X(sel);
		Y = Y(sel);
		Z = round(Z(sel));
		
		figure(600+jj)
		subplot(4,4,ii)
		cla
		hold on
		scatter(X,Y,Z,[.7 .7 1],'s','filled');
		plot(xA,yA,'s','MarkerFaceColor','w','Color',[0 0 .7]);
		
		% visual
		nsamples	= size(samples.pthetaV,1);
		sel			= xSV==sub & xWV==word;
		strV			= ['z_V = ' num2str(zV(sel)) ', N_V = ' num2str(NV(sel)) ,...
			', z_\gamma_V = ' num2str(gV(sel)) ', N_\gamma_V = ' num2str(NgV(sel)) ];

		yV			= zV(sel)./NV(sel);
		nV			= NV(sel);
		z			= NaN(nsamples,1);
		for kk = 1:nsamples
			pthetaV			= samples.pthetaV(kk,sub);
			wthetaV			= samples.wthetaV(kk,word);
			gammaV			= squeeze(samples.gammaV(kk,sub,word));
			
			theta			= pthetaV.*wthetaV;
			theta			= gammaV+(1-gammaV).*theta;
			tmp				= binornd(nV,theta)./nV;
			z(kk)			= tmp;
		end
		u		= unique(z);
		Z		= NaN(size(X));
		for ll = 1:numel(u)
			sel =  z==u(ll);
			Z(ll)	= round(sum(sel)*1000);
		end
		Z = Z./max(Z(:));
		Z = Z*300;
		Y = Y(:);
		Z = Z(:);
		sel = isnan(Z) | Z<=0.9;
		sel = ~sel;
		Y = Y(sel);
		Z = round(Z(sel));
		
		figure(600+jj)
		subplot(4,4,ii)
		scatter(repmat(-23,size(Z)),Y,Z,[1 .7 .7],'s','filled');
		plot([min(xA) max(xA)],[yV yV],'-','MarkerFaceColor','w','Color',[.7 0 0]);
		plot(-23,yV,'s','MarkerFaceColor','w','Color',[.7 0 0]);
		
		% audiovisual
		nsamples = size(samples.pthetaA,1);
		sel		= xSAV==sub & xWAV==word;
% 		strA			= ['z_A = ' num2str(zA(sel)') ', N_A = ' num2str(NA(sel)') ,...
% 			', z_\gamma_A = ' num2str(sum(gA(sel))') ', N_\gamma_A = ' num2str(sum(NgA(sel))') ];

		yAV		= zAV(sel)./NAV(sel);
		xAV		= xSNRAV(sel);
		nAV		= NAV(sel);
		nx		= numel(xAV);
		z		= NaN(nsamples,nx);
% 		for kk = 1:nsamples
% 			pthetaA			= samples.pthetaA(kk,sub);
% 			wthetaA			= samples.wthetaA(kk,word);
% 			pomegaA			= samples.pomegaA(kk,sub);
% 			womegaA			= samples.womegaA(kk,word);
% 			gammaA			= squeeze(samples.gammaA(kk,sub,word));
% 			
% 			theta			= logisticfun(xA,pthetaA+wthetaA,pomegaA+womegaA,0.1);
% 			theta			= gammaA+(1-gammaA).*theta;
% 			tmp				= binornd(nA,theta)./nA;
% 			z(kk,:)			= tmp;
% 		end
% 		X		= xA;
% 		u		= unique(z);
% 		[X,Y]	= meshgrid(X,u);
% 		Z		= NaN(size(X));
% 		for kk	= 1:numel(xA)
% 			for ll = 1:numel(u)
% 				sel =  z(:,kk)==u(ll);
% 				Z(ll,kk)	= round(sum(sel)*1000);
% 			end
% 		end
% 		Z = Z./max(Z(:));
% 		Z = Z*500;
% 		X = X(:);
% 		Y = Y(:);
% 		Z = Z(:);
% 		sel = isnan(Z) | Z<=0.9;
% 		sel = ~sel;
% 		X = X(sel);
% 		Y = Y(sel);
% 		Z = round(Z(sel));
		
		figure(600+jj)
		subplot(4,4,ii)
		hold on
% 		scatter(X,Y,Z,[.7 .7 1],'s','filled');
		plot(xAV,yAV,'s','MarkerFaceColor','w','Color',[0 .7 0]);

		% graphics
		axis square
		box off
		ylim([-0.1 1.5])
		xlim([min(xA)-3 max(xA)+1])
		xlabel('SNR (dB)');
		ylabel('P');
		set(gca,'XTick',xA,'YTick',0:0.2:1,'TickDir','out');
		bf_text(0.1,0.8,strV);
		bf_text(0.1,0.9,strA);
		title(T{word});
	end
end
%%
keyboard

function [samples,stats] = genMCMC(zA,NA,gA,NgA,xSA,xWA,xSNRA,zV,NV,gV,NgV,xSV,xWV,zAV,NAV,gAV,NgAV,xSAV,xWAV,xSNRAV,varargin)
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
parameters		= {'pthetaA','wthetaA','pthetaV','wthetaV','pomegaA','womegaA','gammaA','gammaV','gammaAV'};
writemodel;


% Specify data, as a structure
np = max(xSA);
nw = max(xWA);

dataStruct = struct('zA',zA,'NA',NA,'gA',gA,'NgA',NgA,'xA',xSNRA,'pA',xSA,'wA',xWA,'nTotalA',numel(zA),...
	'zV',zV,'NV',NV,'gV',gV,'NgV',NgV,'pV',xSV,'wV',xWV,'nTotalV',numel(zV),...
	'zAV',zAV,'NAV',NAV,'gAV',gAV,'NgAV',NgAV,'xAV',xSNRAV,'pAV',xSAV,'wAV',xWAV,'nTotalAV',numel(zAV),...
	'np',np,'nw',nw);

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
for ii = 1:nChains
	initsStruct(ii).pthetaA				= zeros(np,1);
	initsStruct(ii).wthetaA				= zeros(nw,1);
	initsStruct(ii).pthetaV				= zeros(np,1);
	initsStruct(ii).wthetaV				= zeros(nw,1);
	initsStruct(ii).pomegaA				= repmat(0.5,np,1);  % ~std
	initsStruct(ii).womegaA				= repmat(0.5,nw,1);  % ~std
	initsStruct(ii).gammaA				= repmat(0.1,np,nw);  % ~std
	initsStruct(ii).gammaV				= repmat(0.1,np,nw);  % ~std
	initsStruct(ii).gammaAV				= repmat(0.1,np,nw);  % ~std
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


parameterNames		= fieldnames(samples); % get all parameter names


%% Save samples
if ~isempty(saveName)
	% 	save([saveName 'Mcmc'],'samples');
end

function writemodel
% Placeholder function to generate JAGS model

%% Model

% 		fun = 'ilogit';
% 		par = '(2*log(1/alpha-1))/omega[s[i]] * (x[i]-theta[s[i]])';
%
% mustr		=	['\t\t f[p[i],w[i]]\t<-    ' fun '( ' par ' ) \r\n'];

%%
str = [
	'\tmodel{ \r\n',...
	'\t# AUDITORY Word performance Is Binomially Distributed \r\n',...
	'\tfor (i in 1:nTotalA){ \r\n',...
	'\t\t zA[i] ~ dbin(thetaA[i],NA[i]) \r\n',...
	'\t\t gA[i] ~ dbin(gammaA[pA[i],wA[i]],NgA[i]) \r\n',...
	'\t\t # Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t\t thetaA[i]  <- (1-gammaA[pA[i],wA[i]])*fA[i] +gammaA[pA[i],wA[i]]\r\n',...
	'\t\t fA[i]	<- ilogit( (2*log(1/alpha-1))/(pomegaA[pA[i]]+womegaA[wA[i]]) * (xA[i]-pthetaA[pA[i]]-wthetaA[wA[i]]) )  \r\n',...
	'\t} \r\n',...
	'\t# VISUAL Word performance Is Binomially Distributed \r\n',...
	'\tfor (i in 1:nTotalV){ \r\n',...
	'\t\t zV[i] ~ dbin(thetaV[i],NV[i]) \r\n',...
	'\t\t gV[i] ~ dbin(gammaV[pV[i],wV[i]],NgV[i]) \r\n',...
	'\t\t # Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t\t thetaV[i] <- (1-gammaV[pV[i],wV[i]])*pthetaV[pV[i]]*wthetaV[wV[i]] +gammaV[pV[i],wV[i]]\r\n',...
	'\t} \r\n',...
	'\t# AUDIOVISUAL Word performance Is Binomially Distributed \r\n',...
	'\tfor (i in 1:nTotalAV){ \r\n',...
	'\t\t zAV[i] ~ dbin(thetaAV[i],NAV[i]) \r\n',...
	'\t\t gAV[i] ~ dbin(gammaAV[pAV[i],wAV[i]],NgAV[i]) \r\n',...
	'\t\t # Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t\t thetaAV[i]  <- (1-gammaAV[pAV[i],wAV[i]])*(fAV[i]+pthetaV[pAV[i]]*wthetaV[wAV[i]]-fAV[i]*pthetaV[pAV[i]]*wthetaV[wAV[i]]) +gammaAV[pAV[i],wAV[i]]\r\n',...
	'\t\t fAV[i]	<- ilogit( (2*log(1/alpha-1))/(pomegaA[pAV[i]]+womegaA[wAV[i]]) * (xAV[i]-pthetaA[pAV[i]]-wthetaA[wAV[i]]) )  \r\n',...
	'\t} \r\n',...
	'\t# Priors For People and Words \r\n',...
	'\tfor (i in 1:np){ \r\n',...
	'\t\t pomegaA[i]\t~ dgamma(1,1)  \r\n',...
	'\t\t pthetaA[i]\t~ dnorm(0,1/(10)^2)\r\n',...
	'\t\t pthetaV[i] ~ dbeta(1,1) \r\n',...
	'\t} \r\n',...
	'\tfor (j in 1:nw){ \r\n',...
	'\t\t womegaA[j]\t~ dgamma(1,1)  \r\n',...
	'\t\t wthetaA[j]\t~ dnorm(0,1/(10)^2)\r\n',...
	'\t\t wthetaV[j] ~ dbeta(1,1) \r\n',...
	'\t} \r\n',...
	'\tfor (i in 1:np){ \r\n',...
	'\t\tfor (j in 1:nw){ \r\n',...
	'\t\t\t gammaA[i,j] ~ dbeta(1,1) \r\n',...
	'\t\t\t gammaV[i,j] ~ dbeta(1,1) \r\n',...
	'\t\t\t gammaAV[i,j] ~ dbeta(1,1) \r\n',...
	'\t\t} \r\n',...
	'\t} \r\n',...
	'\t# constant/default width at F-1(ALPHA) and F-1(1-ALPHA) \r\n',...
	'\t alpha\t\t\t<- 0.1 \r\n',...
	'}\r\n',...
	];
% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);

function [z,N,g,Ng,xSNR,xS,xW,nsnr,usnr] = getdata(M,C,R,T,O,SNR,subjects,words,uwords,mod)

sel			= M==mod;
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
if mod==2
	snr = ones(size(snr));
end
usnr		= unique(snr);
nsnr		= numel(usnr);
z			= NaN(ns,nw,nsnr);
N			= z;
g			= z;
Ng			= z;
xSNR		= z;
xS			= z;
xW			= z;
for ii = 1:ns % for every subject
	for jj = 1:nw % for every word
		for kk = 1:nsnr % for every signal-to-noise ratio
			selsub			= s==us(ii);
			selword			= w==uw(jj);
			selsnr			= snr==usnr(kk);
			sel				= selsub & selword & selsnr';
			
			co				= unique(o(selword)); % current order
			selorder		= o==co;
			selresp			= strcmp(uwords(uw(jj)),r)';
			selstim			= strcmp(uwords(uw(jj)),t)';
			selfalse		= selsub & selresp & ~selstim & selorder' &selsnr';
			
			
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
%%

%%
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

%% 
% [~,~,subs] = unique([xS xW],'rows');
% whos g Ng
% g	= accumarray(subs,g);
% Ng	= accumarray(subs,Ng);
% whos g Ng
function plotthispost(samples)
parameterNames	= fieldnames(samples); % get all parameter names
n = numel(parameterNames);
for parIdx			= 1:n
	a = samples.(parameterNames{parIdx});
	a = a(:);
	subplot(n,n,n*(parIdx-1)+parIdx);
	cla
	plotpost(a);
end

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

function plotparams(samples,T)

figure(100)

% Theta
p			= samples.pthetaA;
w			= samples.wthetaA;
np			= size(p,2);
nw			= size(w,2);

for ii = 1:np
	postSummaryP(ii) = summarizepost(p(:,ii));
end
mu = [postSummaryP.mean];
mP = mu;
E = [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
[~,idx] = sort(mu,2,'descend');
mu		= mu(idx);
E = E(:,idx);
x = 1:np;

subplot(231)
cla
errorpatch(x,-mu,-E);
hold on
plot(x,-mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out','XTick',1:3:np,'XTickLabel',idx(1:3:np),'YTick',-5:5:15,'YTickLabel',5:-5:-15);
box off
xlabel('Subject');
ylabel('Subjective auditory threshold (dB)');
xlim([0 np+1]);
ylim([-5 15])
% horline([0.1  1],'k:');
bf_text(0.05,0.95,'A');
% plot([1 np],[0.1 1],'k:');

for ii = 1:nw
	postSummaryW(ii) = summarizepost(w(:,ii));
end

mu = [postSummaryW.mean];
mW = mu;
E = [[postSummaryW.hdiLow]; [postSummaryW.hdiHigh]];
[~,idx] = sort(mu,2,'descend');
mu		= mu(idx);
E = E(:,idx);
subplot(232)
cla
errorpatch(1:nw,-mu,-E);
hold on
plot(1:nw,-mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out')
words = unique(T);
words = words(idx);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90,'YTick',-5:5:15,'YTickLabel',5:-5:-15);

box off
xlabel('Word');
ylabel('Word auditory threshold (dB)');
xlim([0 nw+1]);
ylim([-5 15])
bf_text(0.05,0.95,'B');

subplot(233)
[mW,mP] = meshgrid(mW,mP);
mu=mP+mW;
mu = -mu;
whos mu
[~,idx1] = sort(mean(mu));
mu		= mu(:,idx1);
[~,idx2] = sort(mean(mu,2));
mu		= mu(idx2,:);
imagesc(mu);
hold on
contour(mu,5,'k');
axis square;
set(gca,'TickDir','out','YDir','normal');
box off
words = unique(T);
words = words(idx1);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);
set(gca,'YTick',1:3:np,'YTickLabel',idx2(1:3:np));


% % Omega
% p			= samples.pomegaA;
% w			= samples.womegaA;
% np			= size(p,2);
% nw			= size(w,2);
%
% for ii = 1:np
% 	postSummaryP(ii) = summarizepost(p(:,ii));
% end
% mu = [postSummaryP.mean];
% E = [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
% [~,idx] = sort(mu);
% mu		= mu(idx);
% E = E(:,idx);
% x = 1:np;
%
% subplot(234)
% cla
% errorpatch(x,mu,E);
% hold on
% plot(x,mu,'ko','MarkerFaceColor','w');
% axis square;
% set(gca,'TickDir','out','XTick',1:3:np,'XTickLabel',idx(1:3:np));
% box off
% xlabel('Subject');
% ylabel('Subjective auditory width (dB)');
% xlim([0 np+1]);
% ylim([-2 10])
% % horline([0.1  1],'k:');
% bf_text(0.05,0.95,'C');
% % plot([1 np],[0.1 1],'k:');
%
% for ii = 1:nw
% 	postSummaryW(ii) = summarizepost(w(:,ii));
% end
%
% mu = [postSummaryW.mean];
% E = [[postSummaryW.hdiLow]; [postSummaryW.hdiHigh]];
% [~,idx] = sort(mu);
% mu		= mu(idx);
% E = E(:,idx);
% subplot(235)
% cla
% errorpatch(1:nw,mu,E);
% hold on
% plot(1:nw,mu,'ko','MarkerFaceColor','w');
% axis square;
% set(gca,'TickDir','out')
% words = unique(T);
% words = words(idx);
% idx = [1 5:5:nw];
% set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);
%
% box off
% xlabel('Word');
% ylabel('Word auditory width (dB)');
% xlim([0 nw+1]);
% ylim([-2 10])
% % horline([0.1 1],'k:');
% % plot([1 50],[0.1 1],'k:');
% bf_text(0.05,0.95,'D');
%
%
% % Gamma
% g			= samples.gammaA;
%
% [~,np,nw]			= size(g);
%
% mu	= NaN(np,nw);
% L	= NaN(np,nw);
% U	= NaN(np,nw);
% for ii = 1:np
% 	for jj = 1:nw
% 		postSummaryG(ii,jj) = summarizepost(squeeze(g(:,ii,jj)));
% 		mu(ii,jj) = postSummaryG(ii,jj).mean;
% 		L(ii,jj) = postSummaryG(ii,jj).hdiLow;
% 		U(ii,jj) = postSummaryG(ii,jj).hdiHigh;
%
% 	end
% end
%
%
% [~,idx1] = sort(mean(mu));
% mu		= mu(:,idx1);
% [~,idx2] = sort(mean(mu,2));
% mu		= mu(idx2,:);
% subplot(233)
% imagesc(mu);
% axis square;
% colorbar;
% caxis([0 0.2])
% set(gca,'TickDir','out','YDir','normal');
% box off
% words = unique(T);
% words = words(idx1);
% idx = [1 5:5:nw];
% set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);
% set(gca,'YTick',1:3:np,'YTickLabel',idx2(1:3:np));


%V

p			= samples.pthetaV;
w			= samples.wthetaV;
np			= size(p,2);
nw			= size(w,2);

for ii = 1:np
	postSummaryP(ii) = summarizepost(p(:,ii));
end
mu = [postSummaryP.mean];
mP = mu;
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
mW = mu;
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

subplot(236)
[mW,mP] = meshgrid(mW,mP);
mu=mP.*mW;
whos mu
[~,idx1] = sort(mean(mu));
mu		= mu(:,idx1);
[~,idx2] = sort(mean(mu,2));
mu		= mu(idx2,:);
imagesc(mu);
hold on
contour(mu,0:0.2:1,'k');
axis square;
% colorbar;
caxis([0 1])
set(gca,'TickDir','out','YDir','normal');
box off
words = unique(T);
words = words(idx1);
idx = [1 5:5:nw];
set(gca,'XTick',idx,'XTickLabel',words(idx),'XTickLabelRotation',90);
set(gca,'YTick',1:3:np,'YTickLabel',idx2(1:3:np));
