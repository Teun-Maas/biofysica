%% Initialization
close all
clearvars

% 	  mu[i] <- alpha*exp(-(zx[i]-zgmu)^2/(2*zgsigma^2))
fun		= @(beta,x)( beta(1).*exp( -1/2 .* (x-beta(2)).^2 ./ beta(3)^2 ) );
x		= linspace(-20,10,5);
nrep	= 2;
x		= repmat(x,1,nrep);

alpha	= 9.3;
theta	= -11;
gsigma	= 6;
sigma	= 0.3;

beta	= [alpha theta gsigma];
y		= fun(beta,x)+sigma*randn(size(x));



scores		= 1:10;
nLevel		= max(scores);
threshold	= 1.5:(nLevel-0.5);

score		= y;
nthresh = numel(threshold);
sel = y<=threshold(1);
N(1) = sum(sel);
score(sel) = 1;
for ii = 2:nthresh
	sel = y<=threshold(ii) & y>threshold(ii-1);
	N(ii) = sum(sel);
	score(sel) = ii;

end
sel = y>threshold(nthresh);
N(nLevel) = sum(sel);
score(sel) = nLevel;


	%%
figure(1)
clf;
plot(x,y,'.','MarkerSize',15);
hold on
plot(x,score,'ks','MarkerSize',15,'LineWidth',2);
nicegraph;
horline(threshold);
horline(scores,'k:');
set(gca,'YTick',scores);
%%
% numSavedSteps	= keyval('numSavedSteps',varargin,5000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
% thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
% burnInSteps		= keyval('burnInSteps',varargin,1000);

beta0		= [alpha theta gsigma];
beta1		= nlinfit(x,score,fun,beta0);
xpredls		= linspace(min(x),max(x),100);
ypredls		= fun(beta1,xpredls);


figure(2)
samples = singleomegafit(x,score,'showDiag',true,'nLevels',nLevel,'numSavedSteps',5000,'burnInSteps',1000);


figure(1)
beta	= [mean(samples.alpha) mean(samples.theta) mean(samples.omega)];
xpred = linspace(min(x),max(x),100);
ypred		= fun(beta,xpred);
hold on
plot(xpred,ypred,'r-','LineWidth',2);

plot(xpredls,ypredls,'k-','LineWidth',2);
%%

%%
thresh		= samples.thresh;
mu			= mean(thresh,2);

figure(142)
clf;
subplot(211);
plot(thresh,mu,'.');
verline(mean(thresh));
xlabel('threshold');

subplot(2,4,5)
plotpost(samples.theta,'xlim',[min(x) max(x)]);
verline(theta);
xlabel('\mu gauss');
title(theta);

subplot(2,4,6)
plotpost(samples.omega,'xlim',[0 max(abs(x))]);
verline(gsigma);
xlabel('\sigma gauss');
title(gsigma);

subplot(2,4,7)
plotpost(samples.alpha,'xlim',[0 nLevel]);
verline(alpha);
xlabel('\alpha gauss');
title(alpha);

subplot(2,4,8)
plotpost(samples.sigma);
verline(sigma);
xlabel('\sigma response');
title(sigma);
