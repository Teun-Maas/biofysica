function SPIN_probsum_model(varargin)
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
% plot_a_vs_v;
plot_av;

function plot_av

datadir			= '/Volumes/mbaudit4/Marc van Wanrooij/data/words';
cd(datadir);
load('spinauditoryMCMC4');
samplesA = samples;
samplesA
ncat = size(samplesA.theta,2)


load('spinvisualMCMC4');
samplesV = samples;


load('spinavMCMC4');
samplesAV = samples;

load('spin');
words = unique(T);



%% Auditory
x		= [-21 -16 -13 -10 -5];
nx		= numel(x);
Apred	= NaN(ncat,nx);
for ii = 1:ncat
	theta = samplesA.theta(:,ii);
	omega = samplesA.omega(:,ii);
	gamma = samplesA.gamma(:,ii);
	n		= numel(theta);
	ypred	= NaN(n,nx);
	for jj = 1:length(gamma)
		ypred(jj,:) = psifun(x,theta(jj),omega(jj),gamma(jj),0,0.1,'function',@logisticfun);
	end
	Apred(ii,:) = mean(ypred);
end

%% Audiovisual
x		= [-21 -16 -13 -10 -5];
nx		= numel(x);
AVpred	= NaN(ncat,nx);
for ii = 1:ncat
	theta = samplesAV.theta(:,ii);
	omega = samplesAV.omega(:,ii);
	gamma = samplesAV.gamma(:,ii);
	n		= numel(theta);
	ypred	= NaN(n,nx);
	for jj = 1:length(gamma)
		ypred(jj,:) = psifun(x,theta(jj),omega(jj),gamma(jj),0,0.1,'function',@logisticfun);
	end
	AVpred(ii,:) = mean(ypred);
end


%% visual
Vpred	= NaN(ncat,1);
for ii = 1:ncat
	ypred = samplesV.theta(:,ii);
	Vpred(ii,:) = mean(ypred);
end
Vpred = repmat(Vpred,1,nx);
%%

ProbSum = Apred+Vpred-Apred.*Vpred;

%%

%%
figure(100)
clf
subplot(231)
plot(Apred,AVpred,'o')
axis square
ylim([-0.1 1.1]);
xlim([-0.1 1.1]);
unityline;
xlabel('P(Correct) A');
ylabel('P(Correct) AV');

subplot(232)
plot(Vpred,AVpred,'o')
axis square
ylim([-0.1 1.1]);
xlim([-0.1 1.1]);
unityline;
xlabel('P(Correct) V');
ylabel('P(Correct) AV');

best = max(Apred,Vpred);
subplot(233)
plot(best,AVpred,'o')
axis square
ylim([-0.1 1.1]);
xlim([-0.1 1.1]);
unityline;
xlabel('P(Correct) Best');
ylabel('P(Correct) AV');

subplot(234)
plot(best,AVpred-best,'o')
axis square
ylim([-0.3 0.3]);
xlim([-0.1 1.1]);
lsline
ylim([-0.3 0.3]);
xlim([-0.1 1.1]);
horline(0,'k:')
xlabel('P(Correct) Best');
ylabel('P(Correct) AV - P(correct) Best');

subplot(235)
plot(ProbSum,AVpred,'o')
hold on
plot(mean(ProbSum,2),mean(AVpred,2),'ko')
xlabel('P(Correct) Summation');
ylabel('P(Correct) AV');

axis square
ylim([-0.1 1.1]);
xlim([-0.1 1.1]);
unityline;




