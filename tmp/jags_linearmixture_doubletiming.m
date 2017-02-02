function jags_linearmixture_doubletiming

%% Initialization
close all


%% The data specification:

delay	= [0 5 10 20 40 80 120 160 320]';
ndelays = numel(delay);
subjects = 1:9;
nsubj = numel(subjects);
for jj = 1:nsubj
	Omega		= NaN(ndelays,1);
	Omegahdi	= NaN(ndelays,2);
	Pi			= Omega;
	Pihdi		= Omegahdi;
	Tau			= Omega;
	Tauhdi		= Omegahdi;
	Tau1			= Omega;
	Tau1hdi		= Omegahdi;
	Beta			= Omega;
	Betahdi		= Omegahdi;
	
	for ii = 1:ndelays
		subj = subjects(jj);
		dT		= delay(ii);
		[y,x,N,y1,x1,N1] = getdata(subj,dT);
% 		plotdata(x,y,x1,y1); % visually check the data
		
		%% JAGS model
		samples		= jagsmodel(y,x,N,y1,x1,N1); % the model & analysis, for afficionados
		%% Check samples
		% 		checkchaincluster(samples);
		samples		= extractchain(samples);
		
		%% Plot relevant parameters
% 				plotparams(samples);
		% 		drawnow
		
		
		%%
		Omega(ii)		= mean(samples.omega);
		Omegahdi(ii,:)	= hdimcmc(samples.omega);
		
		Pi(ii)			= mean(samples.pi(:,1));
		Pihdi(ii,:)		= hdimcmc(samples.pi(:,1));
		
		Tau(ii)		= mean(sqrt(1./samples.tau));
		Tauhdi(ii,:)	= hdimcmc(sqrt(1./samples.tau));
		
		Tau1(ii)		= mean(sqrt(1./samples.tau1));
		Tau1hdi(ii,:)	= hdimcmc(sqrt(1./samples.tau1));
		
		Beta(ii)		= mean(samples.beta1);
		Betahdi(ii,:)	= hdimcmc(samples.beta1);
	end
	
	d = delay';
	d = 1:numel(delay);
	
	figure(42)
	clf
	subplot(221)
	errorpatch(d,Omega,Omegahdi',[.7 0 0]);
	hold on
	plot(d,Omega,'ko','MarkerFaceColor','w','MarkerEdgeColor',[.7 0 0]);
	axis square;
	box off;
	ylim([-0.1 1.1]);

	xlim([0 10]);
	
	set(gca,'TickDir','out',...
		'XTick',1:3:10,...
		'YTick',0:0.2:1);
	xlabel('\Delta T (ms)');
	ylabel('Weight');
	bf_text(0.8,0.8,['Subject: ' num2str(subj)],'HorizontalAlignment','center');
	horline(0.5,'k-');

	
	subplot(222)
	errorpatch(d,Pi,Pihdi',[.7 0 0]);
	hold on
	plot(d,Pi,'ko','MarkerFaceColor','w','MarkerEdgeColor',[.7 0 0]);
	axis square;
	box off;
	ylim([-0.1 1.1]);
 
	xlim([0 10]);
	set(gca,'TickDir','out',...
		'XTick',1:3:10,...
		'YTick',0:0.2:1);
	xlabel('\Delta T (ms)');
	ylabel('P(Response=Weighted)');
	horline(0.5,'k-');

	
	subplot(223)
	errorpatch(d,Beta,Betahdi',[.7 0 0]);
	hold on
	plot(d,Beta,'ko','MarkerFaceColor','w','MarkerEdgeColor',[.7 0 0]);
	axis square;
	box off;
	ylim([-0.1 1.1]);
 
	xlim([0 10]);
	set(gca,'TickDir','out',...
		'XTick',1:3:10,...
		'YTick',0:0.2:1);
	xlabel('\Delta T (ms)');
	ylabel('Gain');
	bf_text(0.8,0.8,['Subject: ' num2str(subj)],'HorizontalAlignment','center');
	horline(0.5,'k-');

	
	subplot(224)
	errorpatch(d,Tau,Tauhdi',[.7 0 0]);
	hold on
	plot(d,Tau,'ko','MarkerFaceColor','w','MarkerEdgeColor',[.7 0 0]);
	errorpatch(d,Tau1,Tau1hdi',[.7 .7 .7]);
	hold on
	plot(d,Tau1,'ko','MarkerFaceColor','w','MarkerEdgeColor',[.7 .7 .7]);
	
	axis square;
	box off;
	ylim([0 15]);
 
	xlim([0 10]);
	set(gca,'TickDir','out',...
		'XTick',1:3:10,...
		'YTick',0:2:14);
	xlabel('\Delta T (ms)');
	ylabel('P(Response=Weighted)');
	horline(0.5,'k-');
	
	
	drawnow
	savegraph([mfilename num2str(subjects(jj))],'eps');
	
	M(jj).omega		= Omega;
	M2(jj).omegahdi	= Omegahdi;
	M(jj).pi		= Pi;
	M2(jj).pihdi		= Pihdi;
	
end
%%
O		= mean([M.omega],2);
hdi		= [M2.omegahdi];
Ohdi	= [mean(hdi(:,1:2:end),2) mean(hdi(:,2:2:end),2)];

P		= mean([M.pi],2);
hdi		= [M2.pihdi];
Phdi	= [mean(hdi(:,1:2:end),2) mean(hdi(:,2:2:end),2)];


figure(99)
clf
subplot(121)
d = delay';
d = 1:numel(d);
errorpatch(d,O',Ohdi')
axis square;
box off;
ylim([-0.1 1.1]);
xlim([-330 330]);
	xlim([0 10]);
set(gca,'TickDir','out',...
	'XTick',1:3:10,...
	'YTick',0:0.2:1);
xlabel('\Delta T (ms)');
ylabel('Weight');
bf_text(0.8,0.8,'Average','HorizontalAlignment','center');
horline(0.5,'k-');
verline(0,'k-');

subplot(122)
errorpatch(d',P',Phdi')
axis square;
box off;
ylim([-0.1 1.1]);
xlim([-330 330]);
	xlim([0 10]);

set(gca,'TickDir','out',...
	'XTick',1:3:10,...
	'YTick',0:0.2:1);
xlabel('\Delta T (ms)');
ylabel('P(Response=Weighted)');
horline(0.5,'k-');
verline(0,'k-');

drawnow
savegraph([mfilename 'average'],'eps');

%%
% for ii = 1:3
% 	figure(ii)
% 	savegraph([mfilename num2str(ii)],'eps');
% end

%%
keyboard
%%

function [y,x,s,N,y1,x1,s1,N1] = getdata(dT)
load('/Users/marcw/DATA/Rachel Grosshardt/Double/A9DataMatrix.mat');
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

function [samples,nChains] = jagsmodel(y,x,N,y1,x1,N1)
%
%
chain_globals;
dataStruct = struct('y',y,...
	'x',x,...
	'N',N,...
	'y1',y1,...
	'x1',x1,...
	'N1',N1);
runjagsMethod = runjagsMethodDefault;
if strcmp(runjagsMethod,'parallel');
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
modelname	= which('linearmixmodel2.txt');
nChains		= 3;
burnInSteps = 1000;
nIter		= 1000;
thinSteps	= 1;
dic			= 0;
parameters	= {'C','z','omega','tau','pi','beta1','tau1'};
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

function plotdata(x,y,x1,y1)

b = regstats(y,x,'linear',{'beta','rsquare'});
b1 = regstats(y1,x1,'linear',{'beta','rsquare'});

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


%% plot
figure(2)
plot(x1,y1,'ko');
str = ['g_1 = ' num2str(b1.beta(2),'%.2f')];
text(30,-50,str,'HorizontalAlignment','center','FontSize',15)
axis square;
axis([-60 60 -60 60]);
xlabel('stimulus x_1 (deg)');
ylabel('response (deg)');
title('Single speaker');
drawnow

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
