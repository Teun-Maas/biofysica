function visualspeechrecognition_singlepar(varargin)
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
datadir			= '/Users/marcw/DATA/Luuk van de Rijt/words';
cd(datadir);
load('spin');
% spin T C SNR L S
roi = 1; % 1 = word, 2 = sentence, 3 = list, 4 = subject

load('spin');
[uwords,~,words]	= unique(T);
[~,~,subjects]		= unique(S); %#ok<*ASGLU>
[ulist,~,list]	= unique(L);
[usent,~,sent]	= unique(Z);

sel			= M==2; % modality 1 = A, 2 = V, 3 = AV
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

us			= unique(s);
ns			= numel(us);
w			= words(sel);
uw			= unique(w);
nw			= numel(uw);

z = NaN(ns,1);
N = z;

for ii = 1:ns
		sel		= s==us(ii);
		z(ii)	= sum(y(sel));
		N(ii)	= sum(sel);
		
end


[~,idx]	= sort(mean(z./N,2));
zsorted		= z(idx,:);
Nsorted		= N(idx,:);

%%
figure(1)
clf
plot(zsorted./Nsorted,'ko-','MarkerFaceColor','w')
axis square
ylim([0 1]);



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
	pause
	
end

%% Extract chain values:
samples = extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

	save(['spinvisualMCMC'  num2str(roi)],'samples');

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



%% Separable rates for Person/subject and Words
figure(777)

% Theta
p			= samples.p;
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
ylim([-0.1 1.1])

% horline([0.1  1],'k:');
bf_text(0.05,0.95,'A');
% plot([1 np],[0.1 1],'k:');

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
% fun				= keyval('function',varargin);
runjagsMethod	= keyval('runjagsMethod',varargin);
dic				= keyval('dic',varargin);
modelname		= fullfile(pwd, 'model.txt');

%% Write the model
% first check parameters to be monitored
parameters		= {'theta','p'};
writemodel;


% Specify data, as a structure
n = length(z);
dataStruct = struct('z',z,'N',N,'n',n);

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
for ii = 1:nChains
	% 	initsStruct(ii).theta			= zeros(np,nw);
	initsStruct(ii).p				= repmat(0.5,n,1);  % ~std
end

%% RUN THE CHAINS
% adaptSteps		= 500;			% Number of steps to 'tune' the samplers.
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethod,'parallel')
	doparallel		= 1; % do use parallelization
	
else
	doparallel		= 0; % do not use parallelization
end


%%
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

%%
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
	'\tfor (i in 1:n){ \r\n',...
	'\t    z[i] ~ dbin(theta[i],N[i]) \r\n',...
	'\t# Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t    theta[i] <- (1-gamma)*p[i] +gamma\r\n',...
	'\t  } \r\n',...
	'\t# Priors For People and Words \r\n',...
	'\tfor (i in 1:n){ \r\n',...
	'\t  p[i] ~ dbeta(1,1) \r\n',...
	'\t} \r\n',...
	'\t  gamma <- 0.1 \r\n',...
	'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);


