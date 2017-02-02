function compare_aud_vs_vis(varargin)
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

load('spinvisualMCMC4');
samplesV = samples;


load('spinavMCMC4');
samplesAV = samples;

load('spin');
words = unique(T);


thetaA = mean(samplesA.theta);
thetaV = mean(samplesV.p);
gammaAV = mean(samplesAV.gamma);
thetaAV = mean(samplesAV.theta);

%%
figure(1)
clf
subplot(321)
h = plot(-thetaA,thetaV,'ko');
hold on
% text(-thetaA,thetaV,num2str([1:18]'),...
% 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% text(thetaV,gammaAV,words,...
% 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% delete(h)
axis square;
box off
xlabel('Threshold A (dB)')
ylabel('P(correct) V')
set(gca,'XTick',10:1:15)
set(gca,'XTickLabel',-10:-1:-15)
set(gca,'FontSize',15)
xlim([9 16]);
ylim([-0.1 1.1]);
unityline;
y = thetaV;
x = -thetaA;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);
regline(beta','k:');
bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1),'%.2f') '-' num2str(pearsonhdi(2),'%.2f') ')']},'FontSize',15);

%%

subplot(323)
h = plot(thetaV,gammaAV,'ko');
hold on
% text(thetaV,gammaAV,num2str([1:18]'),...
% 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% % text(thetaV,gammaAV,words,...
% % 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% delete(h)
axis square;
box off
ylabel('Baseline P(Correct) AV')
xlabel('P(correct) V')
set(gca,'FontSize',15);
xlim([-0.1 1.1]);
ylim([-0.1 1.1]);
unityline;
y = gammaAV;
x = thetaV;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);
regline(beta','k:');
bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1),'%.2f') '-' num2str(pearsonhdi(2),'%.2f') ')']},'FontSize',15);


subplot(324)
h=plot(-thetaA,-thetaAV,'ko');
hold on
% text(-thetaA,-thetaAV,num2str([1:18]'),...
% 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% % text(thetaV,gammaAV,words,...
% % 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% delete(h)
axis square;
box off
xlabel('threshold A (dB)')
ylabel('threshold AV (dB)')
set(gca,'XTick',10:1:15,'XTickLabel',-10:-1:-15,...
	'YTick',10:1:15,'YTickLabel',-10:-1:-15,...
	'FontSize',15);
xlim([9 16])
ylim([9 16])
unityline;
y = -thetaAV;
x = -thetaA;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);
regline(beta','k:');
bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1),'%.2f') '-' num2str(pearsonhdi(2),'%.2f') ')']},'FontSize',15);



subplot(325)
h = plot(thetaV,-thetaAV,'ko');
hold on
% text(thetaV,-thetaAV,num2str([1:18]'),...
% 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% % text(thetaV,gammaAV,words,...
% % 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% delete(h)
axis square
box off
ylabel('threshold AV (dB)')
xlabel('P(correct) V')
set(gca,'YTick',10:1:15,'YTickLabel',-10:-1:-15,...
	'FontSize',15);
ylim([9 16])
xlim([-0.1 1.1]);
horline(mean(-thetaAV),'k:');
verline(mean(thetaV),'k:');
y = -thetaAV;
x = thetaV;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);

regline(beta','k:');

bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1),'%.2f') '-' num2str(pearsonhdi(2),'%.2f') ')']},'FontSize',15);



subplot(326)
h = plot(-thetaA,gammaAV,'ko');
hold on
% text(-thetaA,gammaAV,num2str([1:18]'),...
% 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% % text(thetaV,gammaAV,words,...
% % 	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
% delete(h)
axis square
box off
xlabel('threshold A (dB)')
ylabel('Baseline P(correct) AV')
set(gca,'XTick',10:1:15,'XTickLabel',-10:-1:-15,...
	'FontSize',15);
xlim([9 16])
ylim([-0.1 1.1]);
horline(mean(gammaAV),'k:');
verline(mean(-thetaA),'k:');

y = gammaAV;
x = -thetaA;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);

regline(beta','k:');

bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1),'%.2f') '-' num2str(pearsonhdi(2),'%.2f') ')']},'FontSize',15);

return
% text(-thetaA,thetaV,num2str([1:50]'),'HorizontalAlignment','center','VerticalAlignment','middle');
text(-thetaA,thetaV,words,...
	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
delete(h)
axis square
box off
xlabel('threshold A (dB)')
ylabel('P(correct) V')
set(gca,'XTick',10:1:15,'XTickLabel',-10:-1:-15,...
	'FontSize',15);
xlim([4 21])
ylim([-0.1 1.1]);
horline(mean(thetaV),'k:');
verline(mean(-thetaA),'k:');

%%
y = thetaV;
x = -thetaA;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);

regline(beta','k:');

bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1),'%.2f') '-' num2str(pearsonhdi(2),'%.2f') ')']},'FontSize',15);

savegraph(mfilename,'eps');

function plot_a_vs_v

datadir			= '/Volumes/mbaudit4/Marc van Wanrooij/data/words';
cd(datadir);
load('spinauditoryMCMC1');
samplesA = samples;

load('spinvisualMCMC1');
samplesV = samples;

load('spin');
words = unique(T);


thetaA = mean(samplesA.theta);
thetaV = mean(samplesV.p);


figure(1)
clf
h = plot(-thetaA,thetaV,'ko');
hold on
% text(-thetaA,thetaV,num2str([1:50]'),'HorizontalAlignment','center','VerticalAlignment','middle');
text(-thetaA,thetaV,words,...
	'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
delete(h)
axis square
box off
xlabel('threshold A (dB)')
ylabel('P(correct) V')
set(gca,'XTick',10:1:15,'XTickLabel',-10:-1:-15,...
	'FontSize',15);
xlim([4 21])
ylim([-0.1 1.1]);
horline(mean(thetaV),'k:');
verline(mean(-thetaA),'k:');

%%
y = thetaV;
x = -thetaA;
samples = regjags(y',x');
beta = [mean(samples.beta0) mean(samples.beta1)];
pearsonhdi = hdimcmc(samples.pearson.r);
pearsonmu = mean(samples.pearson.r);

regline(beta','k:');

bf_text(0.9,0.1,...
	{['r^2 = ' num2str(pearsonmu.^2,'%.2f')];['(' num2str(pearsonhdi(1).^2,'%.2f') '-' num2str(pearsonhdi(2).^2,'%.2f') ')']},'FontSize',15);

savegraph(mfilename,'eps');