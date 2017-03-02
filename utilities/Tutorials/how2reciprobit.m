function pa_how2reciprobit
% PA_HOW2RECIPROBIT

% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

close all
clear all
clc

pa_datadir;
load('reactiontime')


%% Histogram
load('reactiontime.mat'); % load the reaction time data into Matlab workspace
rt = RT.easy; % let's take the reaction times (ms) for the "easy" task
x	= 0:50:1000; % reaction time bins (ms)
N	= hist(rt,x);  % count reactions in bin
N	= 100*N./sum(N); % normalization / percentage

figure(1) % graphics
h	= bar(x,N); % histogram/bar plot
% Making figure nicer to look at
set(h,'FaceColor',[.7 .7 .7]); % setting the color of the bars to gray
box off; % remove the "box" around the graph
axis square;
xlabel('Reaction time (ms)'); % and setting labels is essential for anyone to understand graphs!
ylabel('Probability (%)');
xlim([0 1000]);

% optional, save figure
print('-depsc','-painter',[mfilename '1']);

%% Inverse reaction time
rtinv	= 1./rt; % inverse reaction time / promptness (ms-1)

n		= numel(x); % number of bins in reaction time plot
x		= linspace(1/2000,0.01,n); % promptness bins
N		= hist(rtinv,x);
N		= 100*N./sum(N);
figure(2)
h = bar(x*1000,N);
hold on
set(h,'FaceColor',[.7 .7 .7]);
box off
axis square;
xlabel('Promptness (s^{-1})');
ylabel('Probability (%)');
title('Reciprocal time axis');
set(gca,'YTick',0:5:100,'XTick',0:1:8);
xlim([0 8]);

% Does this look like a Gaussian?
% Let's plot a Gaussian curve with mean and standard deviation from the
% promptness data
mu	= mean(rtinv);
sd	= std(rtinv);
y	= normpdf(x,mu,sd);
y	= y./sum(y)*100;
plot(x*1000,y,'ks-','LineWidth',2,'MarkerFaceColor','w');

% optional, save figure
print('-depsc','-painter',[mfilename '2']);

% And test it with a one-sample kolmogorov-smirnov test
% Since this test compares with the normal distribution, we have to
% normalize the data, which we can do with zscore, which is the same as:
% z = (x-mean(x))/std(x)
[h,p]	= kstest(zscore(rt)); % for reaction time
if h
	str		= ['The null hypothesis that the reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
else
	str		= ['The null hypothesis that the reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
end
disp(str); % display the results in command window
[h,p]	= kstest(zscore(rtinv)); % for inverse reaction time
if h
	str		= ['The null hypothesis that the inverse reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
else
	str		= ['The null hypothesis that the inverse reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
end
disp(str); % display the results in command window

%% Cumulative probability
% Normalized scale and nicer shape
x = sort(1000*rtinv);
n = numel(rtinv); % number of data points
y = 100*(1:n)./n; % cumulative probability for every data point
figure(3)
plot(x,y,'k.');
hold on

% Now, plot it cumulative probability in quantiles
% this is easier to compare between different distributions
p		= [1 2 5 10:10:90 95 98 99]/100; % some arbitrary probabilities
q		= quantile(rtinv,p); % instead of hist, let's use quantiles

h = plot(q*1000,p*100,'ko','LineWidth',2,'MarkerFaceColor','r');
hold on
xlabel('Promptness (s^{-1})');
ylabel('Cumulative probability (%)');
title('Cumulative probability plot');
box off
axis square;
set(gca,'YTick',0:10:100,'XTick',0:1:8);
xlim([0 8]);
legend(h,'Quantiles','Location','SE');

% optional, save figure
print('-depsc','-painter',[mfilename '3']);



%% Probit
figure(4)
% raw data
x = -1./sort((rt)); % multiply by -1 to mirror abscissa
n = numel(rtinv); % number of data points
y = pa_probit((1:n)./n); % cumulative probability for every data point converted to probit scale
plot(x,y,'k.');
hold on

% quantiles
p		= [1 2 5 10:10:90 95 98 99]/100;
probit	= pa_probit(p);
q		= quantile(rt,p);
q		= -1./q;
xtick	= sort(-1./(150+[0 pa_oct2bw(50,-1:5)])); % some arbitrary xticks

plot(q,probit,'ko','Color','k','MarkerFaceColor','r','LineWidth',2);
hold on
set(gca,'XTick',xtick,'XTickLabel',-1./xtick);
xlim([min(xtick) max(xtick)]);
set(gca,'YTick',probit,'YTickLabel',p*100);
ylim([pa_probit(0.1/100) pa_probit(99.9/100)]);
axis square;
box off
xlabel('Reaction time (ms)');
ylabel('Cumulative probability');
title('Probit ordinate');

% this should be a straight line
x = q;
y = probit;
b = regstats(y,x);
h = pa_regline(b.beta,'k-');
set(h,'Color','r','LineWidth',2);


% optional, save figure
print('-depsc','-painter',[mfilename '4']);

