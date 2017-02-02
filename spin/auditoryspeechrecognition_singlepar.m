function auditoryspeechrecognition_singlepar(varargin)
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

roi = 1; % 1 = word, 2 = sentence, 3 = list, 4 = subject

loadFlag		= keyval('load',varargin,false);
sampleFlag		= keyval('sample',varargin,true);
datadir			= '/Volumes/mbaudit4/Marc van Wanrooij/data/words';
cd(datadir);

if ~loadFlag
	load('spin');
	[uwords,~,words]	= unique(T);
	[~,~,subjects]		= unique(S);
	[ulist,~,list]	= unique(L);
	[usent,~,sent]	= unique(Z);
	
	sel			= M==1; % modality 1 = A, 2 = V, 3 = AV
	y			= C(sel);
	r			= R(sel);
	t			= T(sel);
	o			= O(sel);
	
	switch roi
		case 4
		s			= subjects(sel);
		case 1
	s			= words(sel);
		case 3
		s			= list(sel);
		case 2
		s			= sent(sel);
	end
	subjects = subjects(sel);
	
	us			= unique(s);
	ns			= numel(us);
	w			= words(sel);
	uw			= unique(w);
	nw			= numel(uw);
	snr			= SNR(sel);
	usnr		= unique(snr);
	nsnr		= numel(usnr);
	
	
	%% Initialization
	z			= NaN(ns,nsnr);
	N			= z;
	xSNR		= z;
	xS			= z;
	
	%% Loop
	for ii = 1:ns % for every word/subject/sentence/list
		for jj = 1:nsnr % for every signal-to-noise ratio
			selsub		= s==us(ii);
			selsnr		= snr==usnr(jj);
			sel			= selsub & selsnr';
			
			z(ii,jj)	= sum(y(sel));
			N(ii,jj)	= sum(sel);
			xSNR(ii,jj)	= usnr(jj);
			xS(ii,jj)	= ii;
		end
	end
	
	
	
% %%
% 	usubjects			= unique(subjects);
% 	nsubjects			= numel(usubjects);
% 	Q = NaN(nsubjects,ns,nsnr);
% 
% 	for kk = 1:nsubjects
% 		for ii = 1:ns % for every subject
% 		for jj = 1:nsnr % for every signal-to-noise ratio
% 			sels		= subjects==usubjects(kk);
% 			selsub		= s==us(ii);
% 			selsnr		= snr==usnr(jj);
% 			whos sels selsub selsnr
% 			sel			= selsub & selsnr' & sels;
% 			
% 			Q(kk,ii,jj)	= sum(sel);
% 
% 		end
% 		end
% 	end
% 	%%
% 			figure(1)
% clf
% cnt = 0;
% 	for ii = 1:nsubjects
% 		cnt = cnt+1;
% 		a = (squeeze(Q(ii,:,:)))==0;
% 		subplot(4,5,cnt)
% 		imagesc(a)
% 		caxis([0 1])
% 		horline(34,'w-')
% 	end
	
	%%
	%% Plot raw data
	close all
	figure(1)
	rho = z./N;
	for ii = 1:nsnr
		r = squeeze(rho(:,ii));
		
		subplot(2,3,ii)
		bar(r);
		box off
		axis square
		caxis([0 1]);
		set(gca,'YDir','normal','TickDir','out');
		
	end
	subplot(2,3,6)
	
	plot(usnr,rho','o-','MarkerFaceColor','w');
	hold on
	xlim([-25 0]);
	ylim([-0.1 1.1]);
	
	
	z		= z(:);
	N		= N(:);
	xSNR	= xSNR(:);
	xS		= xS(:);
	
	sel		= N~=0;
	z		= z(sel);
	N		= N(sel);
	xSNR	= xSNR(sel);
	xS		= xS(sel);
	
	data.hit			= z;
	data.ntrials		= N;
	data.subject		= xS;
	data.SNR			= xSNR; %#ok<STRNU>
	
	save(['spinauditoryraw' num2str(roi)],'data');
elseif loadFlag
	load('spin');
	
	load(['spinauditoryraw' num2str(roi)]);
	z = data.hit; %#ok<NODEF>
	N = data.ntrials;
	xS = data.subject;
	xSNR = data.SNR;
end


if sampleFlag
	%% MCMC
	% 	chain_globals;
	% 	numSavedSteps	= keyval('numSavedSteps',varargin,1000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
	% 	thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
	% 	burnInSteps		= keyval('burnInSteps',varargin,100);
	% 	saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
	% 	nChains			= keyval('nChains',varargin,nChainsDefault);
	% 	runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
	% 	dic				= keyval('dic',varargin,false);
	%
	% 	% graphic flags
	% 	diagFlag		= keyval('showDiag',varargin,true); % show MCMC diagnostics
	% 	postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
	%
	samples = psifit(xSNR,[z N],xS,'lambda',0);
	% 	%% Actual regression
	% 	samples				= genMCMC(z,N,xS,xW,xSNR,...
	% 		'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	% 		'saveName',saveName,'nChains',nChains,...
	% 		'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	% 		'dic',dic);
	
	% % 	%% Thin
	% % 	if thinSteps>1
	% % 		samples = thinchain(samples,thinSteps);
	% % 	end
	%
	%
	% 	%% MCMC diagnostics
	% 	if diagFlag
	% 		plotdiag(samples);
	% 	end
	
	%% Extract chain values:
	% 	samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D
	
	save(['spinauditoryMCMC' num2str(roi)],'samples');
	
	
	
elseif ~sampleFlag
	
	load('spinauditoryMCMC');
end

%%

%% Separable rates for Person/subject and Words
figure(777)

% Theta
p			= samples.theta;
np			= size(p,2);

for ii = 1:np
	postSummaryP(ii) = summarizepost(p(:,ii));
end
mu		= [postSummaryP.mean];
E		= [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
[~,idx] = sort(mu);
mu		= mu(idx);
E		= E(:,idx);
x		= 1:np;

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
ylim([-20 -5])

% horline([0.1  1],'k:');
bf_text(0.05,0.95,'A');
% plot([1 np],[0.1 1],'k:');


% Omega
p			= samples.omega;
np			= size(p,2);

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
ylim([-2 15])
% horline([0.1  1],'k:');
bf_text(0.05,0.95,'C');
% plot([1 np],[0.1 1],'k:');


%%





