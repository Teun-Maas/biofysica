function jags_linearmixture_example

%% Initialization
close all

%% Change these parameters as desired
% you can vary response variability (sd), ratio of WA vs BI
% responses (ratio), and stimulus dominance for bistable percept (theta)

% Noisy
par.ratio	= 2; % ratio of weighted average responses versus bistable responses
par.theta	= 0.5;
par.sd		= 1;
par.wa      = 0.6;

% % Clean
% par.ratio	= 1; % ratio of weighted average responses versus bistable responses
% par.theta	= 0.5;
% par.sd		= 0;
% par.wa      = 0.5;

%% The data specification:
[y,x,N] = getdata(par); 
plotdata(x,y); % visually check the data

%% JAGS model
samples		= jagsmodel(y,x,N); % the model & analysis, for afficionados

%% Check samples
checkchaincluster(samples);
samples		= extractchain(samples);

%% Plot relevant parameters
figure(3)
clf
subplot(131)
plotpost(samples.beta,'xlim',[-0.1 1.1])
axis square;
xlabel('\beta_1');
ylabel('P');

subplot(132)
plotpost(sqrt(1./samples.tau),'xlim',[0 20])
axis square;
xlabel('\sigma (deg)');
ylabel('P');

subplot(133)
plotpost(samples.pi(:,1),'xlim',[-0.1 1.1])
axis square;
xlabel('P(WA/(WA+BI))');
ylabel('P');

%%
for ii = 1:3
	figure(ii)
	savegraph([mfilename num2str(ii)],'eps');
end

function [y,x,N] = getdata(par)
%% Data
x1		= -45:10:45;
x2		= -45:10:45;
[x1,x2] = meshgrid(x1,x2);
x1		= x1(:);
x2		= x2(:);

yavg	= par.wa*x1+(1-par.wa)*x2;
xavg = [x1 x2]; % two single sounds (implicitly weighted average)

ybi		= x1;
Rbin	= logical(binornd(1,par.theta,size(x1)));
ya		= ybi(~Rbin);
yb		= x2(Rbin);
ybi		= [ya;yb];

x		= [x1 x2];
xa		= x(~Rbin,:); % two single sounds (implicitly weighted average)
xb		= x(Rbin,:); % two single sounds (implicitly weighted average)

xbi = [xa;xb];
y	= [repmat(yavg,par.ratio,1); ybi];
x	= [repmat(xavg,par.ratio,1);xbi];
y	= y+par.sd*randn(size(y));

%% Disparity should be >>sd, otherwise impossible to distinguish/cluster
sel = abs(x(:,1)-x(:,2))>0*par.sd;
x = x(sel,:);
y = y(sel,:);

%% Number of observations/responses
N = numel(y);

function [samples,nChains] = jagsmodel(y,x,N)
%
%
chain_globals;
dataStruct = struct('y',y,...
	'x',x,...
	'N',N);
runjagsMethod = runjagsMethodDefault;
if strcmp(runjagsMethod,'parallel');
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
modelname	= which('linearmixmodel.txt');
nChains		= 3;
burnInSteps = 1000;
nIter		= 1000;
thinSteps	= 1;
dic			= 0;
parameters	= {'C','z','beta','tau','pi'};
initsStruct = struct([]);
for ii		= 1:nChains
	initsStruct(ii).tau		= 10;  % ~std
end
% [samples, stats, structArray] = matjags( ...
samples = matjags( ...
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


function checkchaincluster(samples)

%% Check clustering
figure(2)
clf
subplot(321)
beta = samples.beta;
plot(beta');
ylim([-0.1 1.1])
xlabel('sample #');
ylabel('\beta_1');

subplot(322)
tau = samples.tau;
plot(sqrt(1./tau'));
ylim([0 15])
xlabel('sample #');
ylabel('\sigma');

subplot(323)
p1 = squeeze(samples.pi(:,:,1));
plot(p1');
ylim([-0.1 1.1])
xlabel('sample #');
ylabel('P(weighted-average)');

subplot(324)
p2 = squeeze(samples.pi(:,:,2));
plot(p2');
ylim([-0.1 1.1])
xlabel('sample #');
ylabel('P(bistable)');

%% 
C = NaN(3,size(samples.C,3));
for ii = 1:3
	C(ii,:) = mean(samples.C(ii,:,:)==1);
end

subplot(325)
plot(C')
ylim([-0.1 1.1])
xlabel('response #');
ylabel('P(cluster 1)');
%% 
z = NaN(3,size(samples.C,3));
for ii = 1:3
	z(ii,:) = mean(samples.z(ii,:,:)==1);
end

subplot(326)
plot(z')
ylim([-0.1 1.1])
xlabel('response #');
ylabel('P(x_2)');

function plotdata(x,y)

b = regstats(y,x,'linear',{'beta','rsquare'});

%% plot
figure(1)
subplot(3,3,1)
plot(x(:,1),y,'.');
str = ['g_1 = ' num2str(b.beta(2),'%.2f')];
text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
axis square;
xlabel('stimulus x_1 (deg)');
ylabel('response (deg)');

subplot(3,3,2)
plot(x(:,2),y,'.');
str = ['g_2 = ' num2str(b.beta(3),'%.2f')];
text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
axis square;
xlabel('stimulus x_2 (deg)');
ylabel('response (deg)');



subplot(3,3,3)
plot((x(:,1)+x(:,2))/2,y,'.');
str = ['r^2 = ' num2str(b.rsquare,'%.2f')];
text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
axis square;
xlabel('weighted average (deg)');
ylabel('response (deg)');


subplot(3,3,4)
plot((y-(x(:,1)+x(:,2))/2),'.');
xlabel('response #');
ylabel('response - WA (deg)');

subplot(3,3,5)
d1 = y-x(:,1);
d2 = y-x(:,2);
d = [d1 d2];
mn = min(abs(d),[],2);
plot(mn,'.');
xlabel('response #');
ylabel('min|response - BI| (deg)');


subplot(3,3,6)
plot(y-x(:,1),'.');
xlabel('response #');
ylabel('response - WA (deg)');

subplot(3,3,9)
plot(abs(x(:,1)-x(:,2)));
xlabel('response #');
ylabel('|x_1-x_2| (deg)');

subplot(3,3,7)
plot(x(:,1));
xlabel('response #');
ylabel('x_1 (deg)');

subplot(3,3,8)
plot(x(:,2));
xlabel('response #');
ylabel('x_2 (deg)');

drawnow