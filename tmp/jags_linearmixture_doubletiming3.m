function jags_linearmixture_doubletiming

%% Initialization
close all


%% The data specification:

delay	= [0 5 10 20 40 80 120 160 320]';
ndelays = numel(delay);

% [y,x,s,N,y1,x1,s1,N1] = getdata(20);
% 
% d = round(abs(x(:,2)-x(:,1))/5)*5;
% 
% sel = s==9 & d>35;
% Y = y(sel)./d(sel);
% hist(Y,-2:0.2:2);
% 
% return
for ii = 1:ndelays
% for ii = 9;
	ii
	dT		= delay(ii);
	[y,x,s,N,y1,x1,s1,N1] = getdata(dT);
	
% 	sel = ismember(s,4:9);
% 	sel1 = ismember(s1,4:9);
% 	x = x(sel,:);
% 	y = y(sel);
% 	s = s(sel);
% 	N = sum(sel);
% 	
% 	x1 = x1(sel1,:);
% 	y1 = y1(sel1);
% 	s1 = s1(sel1);
% 	N1 = sum(sel1);

% 	plotdata(x,y,s,x1,y1,s1); % visually check the data
	%% JAGS model
	samples		= jagsmodel(y,x,s,N,y1,x1,s1,N1); % the model & analysis, for afficionados
	
	%% Check samples
	% 		checkchaincluster(samples);
	samples		= extractchain(samples);
	
	%% Plot relevant parameters
	% 				plotparams(samples);
	% 		drawnow
	
	S(ii).samples = samples;
end

%%
close all

parName = {'m'};
d		= 1:ndelays;
mu		= NaN(ndelays,1);
hdi		= NaN(ndelays,2);
for ii	= 1:ndelays
	mu(ii)		= mean(S(ii).samples.(parName{1}));
	hdi(ii,:)	= hdimcmc(S(ii).samples.(parName{1}));
end
figure(42)
subplot(121)
sel = ~isnan(mu);
X = d(sel);
Y = mu(sel)';
E = hdi(sel,:)';
errorpatch(X,Y,E);
hold on
plot(d(sel),mu(sel),'ko','MarkerFaceColor','w')
ylim([-0.1 1.1]);
xlim([0 numel(d)+1])
axis square;
box off
horline([0 0.5 1])
set(gca,'TickDir','out',...
	'XTick',1:ndelays,'XTickLabel',delay);
xlabel('Lead (ms)');
ylabel('Leading Weight')


d		= 1:ndelays;
mu		= NaN(ndelays,1);
hdi		= NaN(ndelays,2);
for ii	= 1:ndelays
	a = S(ii).samples.xi(:,1);
	b = S(ii).samples.xi(:,2);
	c = a./(a+b);
	mu(ii)		= mean(c);
	hdi(ii,:)	= hdimcmc(c);
end

subplot(122)
sel = ~isnan(mu);
X = d(sel);
Y = mu(sel)';
E = hdi(sel,:)';
errorpatch(X,Y,E);
hold on
plot(d(sel),mu(sel),'ko','MarkerFaceColor','w')
ylim([-0.1 1.1]);
xlim([0 numel(d)+1])
axis square;
box off
horline([0 0.5 1])
set(gca,'TickDir','out',...
	'XTick',1:ndelays,'XTickLabel',delay);
xlabel('Lead (ms)');
ylabel('P(Weighted Average)');

savegraph('average_mix','png');
%%

function [y,x,s,N,y1,x1,s1,N1] = getdata(dT)
load('/Users/marcw/Dropbox/Work/manuscripts/Rachel PhD/Papers/Paper 2 A9 Double Sounds/Data/A9DataMatrix.mat');
% 1) Trialnumber
% 2) 1snd location (el)
% 3) 2snd location (el)
% 4) el response
% 5) az resp (az location is always zero so not listed in matrix)
% 6) 1.snd
% 7) 2.snd
% 8) subjID (1-9)
% 9) delay (1:17 (-320:320)),
% 10) 1=single sound trials 2= bzz vs gwn trials 3= bzz&gwn (same location)

delay	= -[-320 -160 -120 -80 -40 -20 -10 -5 0 5 10 20 40 80 120 160 320]';


sel			= ismember(delay(Matrix(:,9)),[dT -dT]) & Matrix(:,10)==2;
y			= Matrix(sel,4);
dt			= delay(Matrix(sel,9));
s			= Matrix(sel,8);
x1			= Matrix(sel,2);
x2			= Matrix(sel,3);
x			= [x1 x2];
sel			= dt==-dT;
x(sel,:)	= [x2(sel) x1(sel)];

sel			= dt>=0;
x			= x(sel,:);
y			= y(sel);
s			= s(sel);
% d			= abs(x(:,2)-x(:,1));
% sel			= d>45;
% x			= x(sel,:);
% y			= y(sel);

% if ismember(subj,[8 9])
% sel			= abs(y)>30;
% x			= x(sel,:);
% y			= y(sel);
% end
N			= numel(y);


sel		= ismember(Matrix(:,10),[1 3]);
y1		= Matrix(sel,4);
x1		= Matrix(sel,2);
N1		= numel(y);
s1		= Matrix(sel,8);

function [samples,nChains] = jagsmodel(y,x,s,N,y1,x1,s1,N1)
%
%

[us,~,s] = unique(s); ns = numel(us);
[us1,~,s1] = unique(s1); ns1 = numel(us1);

chain_globals;
dataStruct = struct('y',y,...
	'x',x,...
	's',s,...
	'N',N,...
	'Nsubj',ns,...
	'y1',y1,...
	'x1',x1,...
	's1',s1,...
	'Nsubj1',ns1,...
	'N1',N1);
runjagsMethod = runjagsMethodDefault;
if strcmp(runjagsMethod,'parallel');
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
modelname	= which('linearmixmodel3.txt');
nChains		= 3;
burnInSteps = 1000;
nIter		= 1000;
thinSteps	= 1;
dic			= 0;
parameters	= {'C','z','omega','beta1','tau1',...
	'beta1mu',...
	'm','kappa',...
	'pi','xi'...
	'tau','mtau'};
initsStruct = struct([]);
for ii		= 1:nChains
	initsStruct(ii).tau		= repmat(10,ns,1);  % ~std
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
figure(3)
clf
subplot(321)
omega = samples.omega;
plot(omega');
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

function plotdata(x,y,s,x1,y1,s1)

us = unique(s);
ns = numel(us)

for ii = 1:ns
	sel = s==us(ii);
	
	b = regstats(y(sel),x(sel,:),'linear',{'beta','rsquare'});
	% b1 = regstats(y1,x1,'linear',{'beta','rsquare'});
	
	%% plot
	figure(ii)
	clf
	subplot(1,3,1)
	plot(x(sel,1),y(sel),'.');
	str = ['g_1 = ' num2str(b.beta(2),'%.2f')];
	text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
	axis square;
	xlabel('stimulus x_1 (deg)');
	ylabel('response (deg)');
	
	subplot(1,3,2)
	plot(x(sel,2),y(sel),'.');
	str = ['g_2 = ' num2str(b.beta(3),'%.2f')];
	text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
	axis square;
	xlabel('stimulus x_2 (deg)');
	ylabel('response (deg)');
	
end
%
% subplot(3,3,3)
% plot((x(:,1)+x(:,2))/2,y,'.');
% str = ['r^2 = ' num2str(b.rsquare,'%.2f')];
% text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
% axis square;
% xlabel('weighted average (deg)');
% ylabel('response (deg)');
%
%
% subplot(3,3,4)
% plot((y-(x(:,1)+x(:,2))/2),'.');
% xlabel('response #');
% ylabel('response - WA (deg)');
%
% subplot(3,3,5)
% d1 = y-x(:,1);
% d2 = y-x(:,2);
% d = [d1 d2];
% mn = min(abs(d),[],2);
% plot(mn,'.');
% xlabel('response #');
% ylabel('min|response - BI| (deg)');
%
%
% subplot(3,3,6)
% plot(y-x(:,1),'.');
% xlabel('response #');
% ylabel('response - WA (deg)');
%
% subplot(3,3,9)
% plot(abs(x(:,1)-x(:,2)));
% xlabel('response #');
% ylabel('|x_1-x_2| (deg)');
%
% subplot(3,3,7)
% plot(x(:,1));
% xlabel('response #');
% ylabel('x_1 (deg)');
%
% subplot(3,3,8)
% plot(x(:,2));
% xlabel('response #');
% ylabel('x_2 (deg)');
%
%
% %% plot
% figure(2)
% plot(x1,y1,'ko');
% str = ['g_1 = ' num2str(b1.beta(2),'%.2f')];
% text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
% axis square;
% axis([-60 60 -60 60]);
% xlabel('stimulus x_1 (deg)');
% ylabel('response (deg)');
% title('Single speaker');
% drawnow

function plotparams(samples)
figure(4)
clf
subplot(231)
plotpost(samples.omega,'xlim',[-0.1 1.1])
axis square;
xlabel('\omega');
ylabel('P');

subplot(232)
plotpost(sqrt(1./samples.tau),'xlim',[0 20])
axis square;
xlabel('\sigma (deg)');
ylabel('P');

subplot(233)
plotpost(samples.pi(:,1),'xlim',[-0.1 1.1])
axis square;
xlabel('P(Response =weighted average)');
ylabel('P');

subplot(234)
plotpost(samples.beta1,'xlim',[-0.1 1.1])
axis square;
xlabel('\beta_1');
ylabel('P');

subplot(235)
plotpost(sqrt(1./samples.tau1),'xlim',[0 20])
axis square;
xlabel('\sigma_1 (deg)');
ylabel('P');
