close all
clearvars;

% http://doingbayesiandataanalysis.blogspot.com/2012/06/mixture-of-normal-distributions.html

%% Generate random data from known parameter values:
% set.seed(47405)
trueM1	= 100;
N1		= 500;
trueM2	= 145; % 145 for first example below; 130 for second example
N2		= 500;
trueSD	= 15;

trueM1	= 0;
trueM2	= 0; % 145 for first example below; 130 for second example
N1		= 200;
N2		= 300;
trueSD1	= 1;
trueSD2	= 10;

effsz	= abs( trueM2 - trueM1 ) / trueSD1;
y1 = randn(N1,1);
y1 = (y1-mean(y1))/std(y1) * trueSD1 + trueM1;
y2 = randn(N2,1);
y2 = (y2-mean(y2))/std(y2) * trueSD2 + trueM2;
y = [y1; y2];
N = length(y);




%% Must have at least one data point with fixed assignment 
% to each cluster, otherwise some clusters will end up empty:
Nclust = 2;
clust = NaN(N,1); 
[~,mn] = min(abs(y)); 
[~,mx] = max(y); 

clust(mn) = 1; % smallest value assigned to cluster 1
clust(mx) = 2; % highest value assigned to cluster 2 
% dataList = list(
%     y = y ,
%     N = N ,
%     Nclust = Nclust ,
%     clust = clust ,
%     onesRepNclust = rep(1,Nclust)
% )
dataStruct		= struct('y',y,'N',N,'Nclust',Nclust,'clust',clust,'onesRepNclust',ones(Nclust,1)); % data structure

modelname = which('gmm.txt');


	doparallel		= 1; % do use parallelization

for ii = 1:3
	initsStruct(ii).muOfClust			= [0 mean(y)]; % standardized mean
		initsStruct(ii).tauOfClust			= [1 1./std(y).^2]; % standardized mean

end

parameters		= {'muOfClust','tauOfClust','pClust','clust'};

chain_globals;

%% MCMC
fprintf( 'Running JAGS...\n' );
[samples, stats] = matjags( ...
	dataStruct, ...                     % Observed data
	modelname, ...    % File that contains model definition
	initsStruct, ...                          % Initial values for latent variables
	'doparallel' , doparallel, ...      % Parallelization flag
	'nchains', 3,...              % Number of MCMC chains
	'nburnin', 500,...              % Number of burnin steps
	'nsamples', 1000, ...           % Number of samples to extract
	'thin', 1, ...                      % Thinning parameter
	'dic',false, ...                       % Do the DIC?
	'monitorparams', parameters, ...     % List of latent variables to monitor
	'savejagsoutput',0, ...          % Save command line output produced by JAGS?
	'verbosity',0, ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
	'cleanup',1);                    % clean up of temporary files?


%%
samples = extractchain(samples);

%%
pC1 = samples.pClust(:,1);
pC2 = samples.pClust(:,2);
mC2 = samples.muOfClust(:,2);
mC1 = samples.muOfClust(:,1);
sC1 = 1./sqrt(samples.tauOfClust(:,1));
sC2 = 1./sqrt(samples.tauOfClust(:,2));

c = mean(samples.clust)-1;

figure(1)
clf
subplot(231)

[f,x] = ksdensity(mC1);
hold on

plot(x,f,'LineWidth',2);
hold on
[f,x] = ksdensity(mC2);
plot(x,f,'LineWidth',2);
nicegraph;
xlabel('mean');
verline([trueM1 trueM2]);
subplot(232)
[f,x] = ksdensity(sC1);
plot(x,f,'-','LineWidth',2);

hold on
[f,x] = ksdensity(sC2);
plot(x,f,'-','LineWidth',2);

nicegraph;
xlabel('standard deviation');
verline([trueSD1 trueSD2]);

subplot(233)
[f,x] = ksdensity(pC1);
plot(x,f,'LineWidth',2);
hold on
[f,x] = ksdensity(pC2);
plot(x,f,'LineWidth',2);
nicegraph;
xlabel('proportion');
verline([N1 N2]./(N1+N2));


subplot(234)
plot(mC1,sC1,'.');
hold on
plot(mC2,sC2,'.');
nicegraph;

subplot(236)
[s,idx] = sort(y);
plot(y(idx),c(idx),'k-','MarkerFaceColor','w','LineWidth',2);
nicegraph;
ylabel('probability cluster 2');

subplot(235)
plotpost(y(:));
