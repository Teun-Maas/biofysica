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

cd('/Users/marcw/DATA/Guus van Bentum');
load('data_2AFC_MF.mat'); % data obtained with function XXX

MF		= [5 120 2000];

%% Convert
z		= data(:,3:5);
n		= repmat(120,size(z));
[nSubj,nCond] = size(z);
c		= repmat([1 2 3],nSubj,1);



%% MCMC
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,20000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,5000);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);




%% Actual statistical analysis
samples				= genMCMC(z,n,c,...
	'numSavedSteps',numSavedSteps,'thinSteps',thinSteps,...
	'saveName',saveName,'nChains',nChains,...
	'burnInSteps',burnInSteps,'runJagsMethod',runjagsMethod,...
	'dic',dic);

%% Extract chain values:
samples				= extractchain(samples); % from multiple -dimension matrix to 1- or 2-D

%% Grand average
% Samples contains hierarchy:
% subject pool mean = omega
% subject pool standard deviation: kappa = mu.*(1-mu)./sd.^2 - 1;
% (kappa+1)*sd^2 = mu*(1-mu)
% sd^2 = mu(1-mu)/(kappa+1)

mu			= samples.omega;
muPrior		= samples.omegaPrior;

clear post
mudata	= mean(z./n);
BF = NaN(nCond,1);
figure(1)
clf
for ii = 1:nCond
	subplot(nCond,1,ii)
% 	post(ii) = plotpost(mu(:,ii),'xlab',['\theta_{' num2str(MF(ii)) ' Hz}'],'showCurve',1,'ROPE',[0.45 0.55]);

post(ii) = plotpost(mu(:,ii),'xlab',['\theta_{' num2str(MF(ii)) ' Hz}'],'showCurve',1);

	xlim([0 1]);
	verline(0.5,'k:');
% 	verline(mudata(:,ii),'k-');
% 	plotpost(muPrior(:,ii),'xlab',['\theta_{' num2str(MF(ii)) ' Hz}'],'showCurve',1,'ROPE',[0.45 0.55]);
	
	BF(ii) = bayesfactor(mu(:,ii),muPrior(:,ii),0.5);
	% Bayes Factor10
	% <1 negative (supports H0)
	% 1-3: barely worth mentioning
	% 3-10: substantial
	%10-30: strong
	% 30-100: very strong
	% >100: decisive
	title(['Bayes Factor_{10} = ' num2str(1/BF(ii),2)])
end
figure(2)
clf
subplot(3,1,1)
plotpost(mu(:,2)-mu(:,1),'xlab','\theta_{120 Hz} - \theta_{5 Hz}','showCurve',1);
xlim([-1 1]);
verline(0,'k:');

subplot(3,1,2)
plotpost(mu(:,3)-mu(:,2),'xlab','\theta_{2000 Hz} - \theta_{120 Hz}','showCurve',1);
xlim([-1 1]);
verline(0,'k:');
subplot(3,1,3)
plotpost(mu(:,3)-mu(:,1),'xlab','\theta_{2000 Hz} - \theta_{120 Hz}','showCurve',1);
xlim([-1 1]);
verline(0,'k:');

%% Individual data
theta	= samples.theta;
mutheta = NaN(nSubj,nCond);


for ii = 1:nCond
	mutheta(:,ii) = mean(squeeze(theta(:,:,ii)));
end

figure(3)
clf
plot(mutheta','k-','Color',[.7 .7 .7])
hold on
errorbar(1:3,[post.mean],[post.mean]-[post.hdiLow],[post.hdiHigh]-[post.mean],...
	'ko-','MarkerSize',10,'MarkerFaceColor','w','LineWidth',2);
axis square
box off
xlim([0.5 3.5])
set(gca,'TickDir','out',...
	'XTick',1:3,'XTickLabel',[5 120 2000],...
	'YTick',0:0.2:1);
ylim([-0.1 1.1]);

ylabel('P(\theta)');
xlabel('Modulation frequency (Hz)');
horline(0.5,'k:');
h = text(1:3,ones(3,1),num2str(1./BF,2),'HorizontalAlignment','center');
% set(h,

%%
print('-depsc','-painters',mfilename); % Mac?
% print('-depsc','-painter',mfilename); % Windows?

%%
keyboard

function [samples,stats] = genMCMC(z,N,c,varargin)
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
	'thetaPrior','omegaPrior','kappaPrior',};
writemodel;

[Nsubj,Ncond]	= size(z);
dataStruct		= struct('z',z,'N',N,'Nsubj',Nsubj,'c',c,'Ncond',Ncond);

%% INTIALIZE THE CHAINS.
initsStruct = struct([]);
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
for ii = 1:nChains
	initsStruct(ii).theta			= repmat(0.5,Nsubj,Ncond);
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

%%

%%
str = [
	'model{ \r\n',...
	'\tfor ( s in 1:Nsubj ) { \r\n',...
	'\tfor ( c in 1:Ncond ) { \r\n',...
	'\t\t\tz[s,c] ~ dbin(theta[s,c], N[s,c] ) \r\n',...
	'\t\t\ttheta[s,c] ~ dbeta( omega[c]*(kappa[c]-2)+1 , (1-omega[c])*(kappa[c]-2)+1 ) \r\n',...
	'\t\t\tthetaPrior[s,c] ~ dbeta( omegaPrior[c]*(kappaPrior[c]-2)+1 , (1-omegaPrior[c])*(kappaPrior[c]-2)+1 ) \r\n',...
	'\t\t} \r\n',...
	'\t} \r\n',...
	'\tfor ( c in 1:Ncond ) { \r\n',...
	'\tomega[c] ~ dbeta(1,1) # broad uniform \r\n',...
	'\tkappa[c] <- kappaMinusTwo[c] + 2\r\n',...
	'\tkappaMinusTwo[c] ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10\r\n',...
	'\t\t} \r\n',...
	'\tfor ( c in 1:Ncond ) { \r\n',...
	'\tomegaPrior[c] ~ dbeta(1,1) # broad uniform \r\n',...
	'\tkappaPrior[c] <- kappaMinusTwoPrior[c] + 2\r\n',...
	'\tkappaMinusTwoPrior[c] ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10\r\n',...
	'\t\t} \r\n',...
	'}\r\n',...
	];


% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);


