function visualspeechrecognition_singlepar(varargin)
% MCMC = BINORATE(Z,N,S)

%% Initialization
close all
cmap		= statcolor(64,[],[],[],'def',6);

% graphic flags
diagFlag		= keyval('showDiag',varargin,false); % show MCMC diagnostics
postFlag		= keyval('showPost',varargin,false); % show posterior estimate distribution
roi				= 1;

%% Load data
[y,subjects,words,uwords]	= loaddata;

%% Arrays
[z,N,p]						= getarray(subjects,words,y);
[zroi,Nroi]					= getroi(roi,subjects,words,y);

%% Separability
% % just a test on separability
% % https://xcorr.net/2011/09/20/using-the-svd-to-estimate-receptive-fields/
% % mup = mean(mu(:));
% % mup=0
% % [u,s,v] = svd(reshape(mu,50,18)-mup);
% 
% mup = mean(p(:));
% mup=0
% [u,s,v] = svd(p-mup);
% lambda = diag(s);
% 
% s1 = zeros(size(z));
% s1(1) = lambda(1);
% % s1(2,2) = lambda(2);
% % s1(3,3) = lambda(3);
% 
% p1 = u*s1*v'+mup;
% 
% close all
% 
% subplot(321)
% imagesc(p)
% caxis([0 1]);
% nicegraph;
% 
% subplot(322)
% plot(-v(:,1)*sqrt(lambda(1)))
% hold on
% plot(mean(p));
% 
% subplot(323)
% plot(-u(:,1)*sqrt(lambda(1)))
% hold on
% plot(mean(p,2));
% 
% subplot(324)
% imagesc(p1);
% caxis([0 1]);
% nicegraph;
% r = corrcoef(p,p1);
% r = r(2)^2;
% title(num2str(r,'%0.2f'))
% 
% subplot(325)
% plot(lambda./sum(lambda),'ko-');
% nicegraph

%% Graph
graph1(p,uwords,zroi,Nroi);


%% Actual regression
samples				= genMCMC(z,N,'model','model2.txt');

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
% Theta
postp			= samples.p;
np				= size(postp,2);
% theta			= samples.p;
for ii = 1:np
	postSummaryP(ii) = summarizepost(postp(:,ii));
end
mu		= [postSummaryP.mean];
E		= [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
[~,idx] = sort(mu);
x		= 1:np;


%%

figure(3)
clf

postp			= samples.p;
np				= size(postp,2);
% theta			= samples.p;
for ii = 1:np
	postSummaryP(ii) = summarizepost(postp(:,ii));
end
mu		= [postSummaryP.mean];
E		= [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
[~,idx] = sort(mu);
x		= 1:np;

subplot(131)
cla
errorpatch(x,mu,E);
hold on
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'XTick',1:3:np,'XTickLabel',uwords(1:3:np),...
	'XTickLabelRotation',45);
box off
xlabel('Word');
ylabel('Visual word recognition rate p');
xlim([0 np+1]);
ylim([-0.1 1.1])
bf_text(0.05,0.95,'A','FontSize',15);
horline([0 1],'k-');
set(get(gca,'Xaxis'),'MinorTickValues',1:50)
nicegraph;


postp			= samples.q;
np				= size(postp,2);
% theta			= samples.p;
clear postSummaryP
for ii = 1:np
	postSummaryP(ii) = summarizepost(postp(:,ii));
end
mu		= [postSummaryP.mean];
E		= [[postSummaryP.hdiLow]; [postSummaryP.hdiHigh]];
[~,idx] = sort(mu);
x		= 1:np;

subplot(132)
cla
errorpatch(x,mu,E);
hold on
plot(x,mu,'ko','MarkerFaceColor','w');
axis square;
set(gca,'TickDir','out','XTick',1:3:np,'XTickLabel',idx(1:3:np));
box off
xlabel('Subject');
ylabel('Subjective visual recognition rate q');
xlim([0 np+1]);
ylim([-0.1 1.1])
bf_text(0.05,0.95,'B','FontSize',15);
horline([0 1],'k-');
nicegraph;
hA=gca;
set(get(gca,'Xaxis'),'MinorTickValues',1:50)
% set(get(gca,'Yaxis'),'MinorTickValues',linspace(0,1,19));

mu			= squeeze(mean(samples.theta));
ax2 = subplot(133);
imagesc(mu);
hold on
% contour(mu,sort([0.05 0:0.25:1 0.95]),'w');
yi = linspace(1,50,8);
xi = linspace(1,18,18);
set(gca,'YDir','normal',...
	'XTick',xi,'YTick',yi);
ylim([0.5 50.5]);
xlim([0.5 18.5]);
axis square;
bf_text(0.95,0.95,'C','FontSize',15,'Color','w');
colormap(ax2,1-gray);
set(gca,'YTick',1:3:50,'YTickLabel',uwords(1:3:50),...
	'XTick',1:np,'XTickLabel',idx(1:np),...
	'XTickLabelRotation',0);
title({['Visual recognition rate'],['(1-\gamma)\cdotp\cdotq+\gamma']});
set(get(gca,'Xaxis'),'MinorTickValues',1:50);
set(get(gca,'Yaxis'),'MinorTickValues',1:50);
xlabel('Subject');
ylabel('Word');
nicegraph;
% colorbar
caxis([0.1 1]);

savegraph('lipreading','eps');
%%
figure(5)
clf

mu			= p;
ax2 = subplot(133);
imagesc(mu)
hold on
% contour(mu,sort([0.05 0:0.25:1 0.95]),'w');
yi =linspace(1,50,8);
xi = linspace(1,18,18);
set(gca,'YDir','normal',...
	'XTick',xi,'YTick',yi);
caxis([0 1]);
ylim([0.5 50.5]);
xlim([0.5 18.5]);
axis square;
nicegraph;
colormap(ax2,cmap);

mu			= squeeze(mean(samples.theta));
est = (z+1)./(N+2);
bubbleplot(est(:),mu(:),'Xwidth',0.1,'Ywidth',0.1);
hold on
plot(est(:),mu(:),'k.');

ylim([-0.1 1.1]);
xlim([-0.1 1.1]);
unityline;
axis square;
r = corrcoef(est,mu);
r = r(2)^2;
str = ['r^2 = ' num2str(r,'%0.2f')];
title(str)
nicegraph;

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
figure(4)
plot(p+0.01*randn(size(p)),0.9*mu+0.1+0.01*randn(size(p)),'k.')
axis square
unityline
nicegraph

%%
n	= mean(N(:));
x	= 0:n;
pmu = mean(postp);
L	= binopdf(z(:),n,pmu(:));

p = log(numel(z)*18);
k = numel(z);
BIC = [sum(-2*log(L))+k*p]
AIC = [sum(-2*log(L))+2*p]
%%

keyboard

function [samples,stats] = genMCMC(z,N,varargin)
% B = GENMCMC(X,Y)
%
% Generate MCMC chains

%% initialization
% most is already defined in main function, so no defaults
chain_globals;
numSavedSteps	= keyval('numSavedSteps',varargin,5000); % number of saved MCMC samples. How many you need depend on autocorrelation (effective sample size>10000), convergence (shrink factor<1.1), etc
thinSteps		= keyval('thinSteps',varargin,1); % 1/proportion MCMC samples thrown away
burnInSteps		= keyval('burnInSteps',varargin,500);
saveName		= keyval('saveName',varargin,'Hier-PsychometricCurve-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);
runjagsMethod	= keyval('runjagsMethod',varargin,runjagsMethodDefault);
dic				= keyval('dic',varargin,false);

mdl				= keyval('model',varargin,'model1.txt');

strcmp(mdl,'model1.txt')
modelname		= fullfile(pwd, mdl);

%% Write the model
% first check parameters to be monitored
% Specify data, as a structure
initsStruct = struct([]);

switch mdl
	case 'model1.txt'
		parameters		= {'theta','p'};
		n = length(z);
		dataStruct = struct('z',z,'N',N,'n',n);
		for ii = 1:nChains
			initsStruct(ii).p				= repmat(0.5,n,1);  % ~std
			initsStruct(ii).q				= repmat(0.5,n,1);  % ~std
		end
	case 'model2.txt'
		[np,nq]= size(z);
		parameters		= {'theta','p','q'};
		dataStruct = struct('z',z,'N',N,'np',np,'nq',nq);
		for ii = 1:nChains
			initsStruct(ii).p				= repmat(0.5,np,1);  % ~std
			initsStruct(ii).q				= repmat(0.5,nq,1);  % ~std
		end
end
writemodel;



%% INTIALIZE THE CHAINS.
% % because covariates/explanatory variables/x data are standardized,
% % means(x) can be set to zero, and stds(x) to 1.
% % Probabilities (guess and lambda) are determined from dependent variables/
% % y data
	% 	initsStruct(ii).theta			= zeros(np,nw);


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
parameterNames		= fieldnames(samples); % get all parameter names


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
fid			= fopen('model1.txt','w');
fprintf(fid,str);
fclose(fid);


%% Model
str = [
	'\tmodel{ \r\n',...
	'\t# Correctness Of Each Answer Is Bernoulli Trial \r\n',...
	'\tfor (i in 1:np){ \r\n',...
	'\tfor (j in 1:nq){ \r\n',...
	'\t    z[i,j] ~ dbin(theta[i,j],N[i,j]) \r\n',...
	'\t# Probability Correct Is Product Of Word By Person Rates \r\n',...
	'\t    theta[i,j] <- (1-gamma)*p[i]*q[j] +gamma\r\n',...
	'\t  } \r\n',...
	'\t  } \r\n',...
	'\t# Priors For People and Words \r\n',...
	'\tfor (i in 1:np){ \r\n',...
	'\t  p[i] ~ dbeta(1,1) \r\n',...
	'\t} \r\n',...
	'\tfor (i in 1:nq){ \r\n',...
	'\t  q[i] ~ dbeta(1,1) \r\n',...
	'\t} \r\n',...
	'\t  gamma <- 0.1 \r\n',...
	'}\r\n',...
	];



% Write the modelString to a file, using Matlab commands:
fid			= fopen('model2.txt','w');
fprintf(fid,str);
fclose(fid);

function [zsorted,Nsorted] = getroi(roi,subjects,words,y)
switch roi
	case 1
		s = words;
	case 2
		s = subjects;
end
us		= unique(s);
ns		= numel(us);

z		= NaN(ns,1);
N		= z;
for ii = 1:ns
	sel		= s==us(ii);
	z(ii)	= sum(y(sel));
	N(ii)	= sum(sel);
	
end
[~,idx]		= sort(mean(z./N,2));
zsorted		= z(idx,:);
Nsorted		= N(idx,:);

function [y,subjects,words,uwords] = loaddata
datadir			= '/Users/marcw/DATA/Luuk van de Rijt/words';
cd(datadir);
load('spin'); % spin T C SNR L S
[uwords,iwords,words]	= unique(T);
[~,~,subjects]		= unique(S); %#ok<*ASGLU>

sel					= M==2; % modality 1 = A, 2 = V, 3 = AV
y					= C(sel);
subjects			= subjects(sel);
words				= words(sel);

function [z,N,p] = getarray(subjects,words,y)

[u,idx,idxb]	= unique([subjects words],'rows');
z				= accumarray(idxb,y,[],@sum);
N				= accumarray(idxb,ones(size(y)),[],@sum);
z				= reshape(z,50,18);
N				= reshape(N,50,18);
p				= z./N;

mu				= mean(p);
[~,idx]			= sort(mu);
z				= z(:,idx);
N				= N(:,idx);
p				= p(:,idx);

mu				= mean(p,2);
[~,idx]			= sort(mu);
z				= z(idx,:);
N				= N(idx,:);
p				= p(idx,:);

function graph1(p,uwords,zroi,Nroi)

cmap		= statcolor(64,[],[],[],'def',6);
cmap		= cmap(:,[3 1 2]);

figure(1)
clf
subplot(131)
plot(p)
hold on
plot(mean(p,2),'ko','MarkerFaceColor','w','MarkerSize',10)
axis square
ylim([0 1]);
xlabel('Word');
ylabel('P(Correct)');
xlim([0 51])
horline(0.1);
nicegraph;
xi			= linspace(1,50,8);
set(gca,'XTick',xi,'XTickLabel',uwords(xi))
title('Lip reading');

subplot(132)
plot(p')
hold on
plot(mean(p),'ko','MarkerFaceColor','w','MarkerSize',10)
axis square
ylim([0 1]);
xlabel('Subject');
ylabel('P(Correct)');
xlim([0 19])
horline(0.1);
nicegraph;
set(gca,'XTick',1:18);
title('Lip reading');

ax2 = subplot(133);
imagesc(p)
yi =linspace(1,50,8);
xi = linspace(1,18,18);
set(gca,'YDir','normal',...
	'XTick',xi,'YTick',yi);
caxis([0 1]);
ylim([0.5 50.5]);
xlim([0.5 18.5]);
axis square;
nicegraph;
colormap(ax2,cmap);


figure(2)
clf
plot(zroi./Nroi,'ko','MarkerFaceColor','w','MarkerSize',10)
axis square
ylim([0 1]);
xlabel('Word');
ylabel('P(Correct)');
xlim([0 51])
horline(0.1);
nicegraph;
xi =linspace(1,50,8);
set(gca,'XTick',xi,'XTickLabel',uwords(xi))
title('Lip reading');