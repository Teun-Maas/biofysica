function ss_gifford_regresspower
% % Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa-Power.R
close all; % This closes all of Matlab's graphics windows.
clc;

seed = 47405; rng(seed); % Optional, merely for replicability.

singleFlag	= false;
diagFlag	=false;
postFlag	= false;
densFlag = false;
% Load the functions genMCMC, smryMCMC, and plotMCMC:
% (This also sources DBDA2E-utilities.R)
% source('Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa.R')

%  Specify idealized hypothesis:
% idealGroupPar0		= [0.3 0.1 0 5]; % ideal group parameters, Cond 0 & Cond 1 Gain Mean & SD, Bias Mean & SD
% idealGroupPar1		= [0.6 0.1 0 5]; % ideal group parameters, Cond 0 & Cond 1 Gain Mean & SD, Bias Mean & SD
% idealNsubj			= 20; % more subjects => higher confidence in hypothesis
% idealNtrl			= 20; % more trials => higher confidence in hypothesis
% idealSigma			= 5; % does not vary per subject
%
% [xVec,cVec,yVec,sVec] = simdata(idealGroupPar0,idealGroupPar1,idealSigma,idealNsubj,idealNtrl);
[xVec,cVec,yVec,sVec] = getdata;

%% Analyse per subject
if singleFlag
	us = unique(sVec);
	ns = numel(us);
	uc = unique(cVec);
	nc = numel(uc);
	for ii = 1:ns
		for jj = 1:nc
			sel = sVec==us(ii) & cVec==uc(jj);
			sum(sel)
			if sum(sel)
				x = xVec(sel);
				y = yVec(sel);
				figure(100+ii)
				subplot(1,2,jj)
				plotloc(x,y);
				drawnow
				
			end
			savegraph(['giffordsubject' num2str(ii)],'png');
		end
	end
end

%% Run Bayesian analysis on idealized data:
mcmc			= hxregjags(xVec,cVec,yVec,sVec);

% keyboard
%% densities
if densFlag
xmcmc = extractchain(mcmc);
	parameterNames	= fieldnames(xmcmc); % get all parameter names
	for parIdx			= 1:numel(parameterNames)
		parName = parameterNames{parIdx};
		s =  xmcmc.(parName);
		nsub = size(s,2);
		if nsub>1
			for ii = 1:nsub
				figure(parIdx)
				subplot(3,5,ii)
		plotpost(s(:,ii));
		title(parName)
			end
		end
	end
end
	
	%% diagnose
if diagFlag
	parameterNames	= fieldnames(mcmc); % get all parameter names
	for parIdx			= 1:numel(parameterNames)
		parName = parameterNames{parIdx};
		s =  mcmc.(parName);
		if ~ismatrix(s)
			[~,~,n] = size(s);
			for sIdx = 1:n
				a = squeeze(s(:,:,sIdx));
				newparName = [parName '_' num2str(sIdx)];
				samples.(newparName) = a;
			end
		else
			samples.(parName) = mcmc.(parName);
			
		end
	end
	
	
	parameterNames	= fieldnames(samples); % get all parameter names
	for parIdx			= 1:numel(parameterNames)
		figure(100+parIdx)
		clf;
		diagmcmc(samples,'parName',parameterNames{parIdx});
		drawnow
			savegraph(['gifforddiag' parameterNames{parIdx}],'png');

	end
	
end



%%

mcmc = extractchain(mcmc);

if postFlag
	%% Bias
	figure(200)
	clf
	m = mcmc.beta1mu;
	s = mcmc.beta1sigma;
	y = m + s.*randn(size(m));
	y = zscore(y);
	y = mean(m)+mean(s)*y;
	subplot(331)
	plotpost(y,'xlim',[-90 90])
	axis square
	
	subplot(332)
	plotpost(m,'xlim',[-90 90])
	axis square
	
	subplot(333)
	plotpost(mcmc.beta1(:),'xlim',[-90 90])
	axis square
	
	
	m = mcmc.beta1mu+mcmc.beta2mu;
	s1 = mcmc.beta1sigma;
	s2 = mcmc.beta2sigma;
	y = m + s1.*randn(size(m))+ s2.*randn(size(m));
	% y = zscore(y);
	% y = mean(m)+mean(s)*y;
	subplot(334)
	plotpost(y,'xlim',[-90 90])
	axis square
	
	subplot(335)
	plotpost(m,'xlim',[-90 90])
	axis square
	
	subplot(336)
	plotpost(mcmc.beta1(:)+mcmc.beta2(:),'xlim',[-90 90])
	axis square
	
	
	m = mcmc.beta2mu;
	s = mcmc.beta2sigma;
	y = m + s.*randn(size(m));
	% y = zscore(y);
	% y = mean(m)+mean(s)*y;
	subplot(337)
	plotpost(y,'xlim',[-90 90])
	axis square
	
	subplot(338)
	plotpost(m,'xlim',[-90 90])
	axis square
	
	subplot(339)
	plotpost(mcmc.beta2(:),'xlim',[-90 90])
	axis square
	savegraph('giffordbiaspost','png');
	
	
	
	%% Gain
	figure(201)
	clf
	m = mcmc.beta3mu;
	s = mcmc.beta3sigma;
	y = m + s.*randn(size(m));
	y = zscore(y);
	y = mean(m)+mean(s)*y;
	subplot(331)
	plotpost(y,'xlim',[-1 2])
	axis square
	
	subplot(332)
	plotpost(m,'xlim',[-1 2])
	axis square
	
	subplot(333)
	plotpost(mcmc.beta3(:),'xlim',[-1 2])
	axis square
	
	
	m = mcmc.beta3mu+mcmc.beta4mu;
	s1 = mcmc.beta3sigma;
	s2 = mcmc.beta4sigma;
	y = m + s1.*randn(size(m))+ s2.*randn(size(m));
	% y = zscore(y);
	% y = mean(m)+mean(s)*y;
	subplot(334)
	plotpost(y,'xlim',[-1 2])
	axis square
	
	subplot(335)
	plotpost(m,'xlim',[-1 2])
	axis square
	
	subplot(336)
	plotpost(mcmc.beta3(:)+mcmc.beta4(:),'xlim',[-1 2])
	axis square
	
	
	m = mcmc.beta4mu;
	s = mcmc.beta4sigma;
	y = m + s.*randn(size(m));
	% y = zscore(y);
	% y = mean(m)+mean(s)*y;
	subplot(337)
	plotpost(y,'xlim',[-1 2])
	axis square
	
	subplot(338)
	plotpost(m,'xlim',[-1 2])
	axis square
	
	subplot(339)
	plotpost(mcmc.beta4(:),'xlim',[-1 2])
	axis square
	savegraph('giffordgainpost','png');
end

%% Simulation
% Specify sample size for each simulated data set:
Nsubj				= 20; Nrep = 1;  % 150 trials total
% Specify the number of simulated experiments:
n					= numel(mcmc.beta1(:));
nSimulatedDataSets	= min([500,n]); % An arbitrary large number.
% Run the simulated experiments:
simCount			= 0;
if exist('goalTally','var')
	clear goalTally % in case previously run from here down
end
% goalTally = struct([]);
if ~exist('goalTally','var') 		   % if goalTally does not exist, create it
	goalTally = [];
end


%%
for simIdx = 1:nSimulatedDataSets
	simCount	= simCount+1;
	disp(['\n\n==================== Simulation ' num2str(simCount) ' of ' num2str(nSimulatedDataSets) '====================\n\n']);
	% Generate random beta1/0mu and beta1/0sd  for group distribution:

	groupPar.beta1mu = mcmc.beta1mu(simIdx);
	groupPar.beta2mu = mcmc.beta2mu(simIdx);
	groupPar.beta3mu = mcmc.beta3mu(simIdx);
	groupPar.beta4mu = mcmc.beta4mu(simIdx);

	groupPar.beta1sigma = mcmc.beta1sigma(simIdx);
	groupPar.beta2sigma = mcmc.beta2sigma(simIdx);
	groupPar.beta3sigma = mcmc.beta3sigma(simIdx);
	groupPar.beta4sigma = mcmc.beta4sigma(simIdx);

	groupPar.sigma = mcmc.sigma(simIdx);
	
	[xVec,cVec,yVec,sVec,Par] = simdata(groupPar,Nsubj,Nrep);
	
	
	%%
	
	%%
	% 	%% Generate random data based on parameter value:
	% 	b1				= normrnd(genBeta1mu,genBeta1sigma,[Nsubj,1]); % gain
	% 	b0				= normrnd(genBeta0mu,genBeta0sigma,[Nsubj,1]); % bias
	%
	%
	% 	xVec			= transpose(round(linspace(-90,90,Ntrl)/15)*15);
	% 	yVec			= NaN(Ntrl,Nsubj);
	% 	for sIdx = 1:Nsubj
	% 		yVec(:,sIdx) = b1(sIdx).*xVec+b0(sIdx)*randn(1);
	% 	end
	% 	sVec			= repmat(1:Nsubj,[Ntrl,1]);
	% 	xVec			= repmat(xVec,[1,Nsubj]);
	%
	% 	yVec			= yVec(:);
	% 	sVec			= sVec(:);
	% 	xVec			= xVec(:);
	
	
	dataStruct.y	= yVec;
	dataStruct.s	= sVec;
	dataStruct.x	= xVec;
	dataStruct.c	= cVec;
	dataStruct.p	= Par;
	% Do Bayesian analysis on simulated data:
	goalAchieved = goalAchievedForSample(dataStruct);
	% Tally the results:
	goalTally = [goalTally goalAchieved];
end

%% For each goal...
% Extract the goal name for subsequent display:
goalName = fieldnames(goalTally);
for goalIdx = 1:numel(goalName)
	goal = [goalTally.(goalName{goalIdx})];
	% Compute number of successes:
	goalHits = sum(goal);
	% Compute number of attempts:
	goalAttempts = numel(goal);
	% Compute proportion of successes:
	goalEst = goalHits/goalAttempts;
	% Compute HDI around proportion:
	goalEstHDI = hdicdf('beta',0.95,1e-08,1+goalHits,1+goalAttempts-goalHits);
	% Display the result:
	str = [goalName{goalIdx} ': Est.Power=' , num2str(round(goalEst,3)) '; Low Bound=' , num2str(round(goalEstHDI(1),3)),...
		'; High Bound=' , num2str(round(goalEstHDI(2),3))];
	disp(str)
end


%% Subfunctions

% Define function that assays goal achievement for a single set of data:
function goalAchieved = goalAchievedForSample(data)
% Generate the MCMC chain:
tic
mcmc		= hxregjags(data.x,data.c,data.y,data.s);
toc
% Convert coda object to matrix for convenience:
% Specify criteria for goals:

goalgain = 0.1;
goalbias = 5;
goalbiasnarrow = 5;
goalgainnarrow = 0.1;
mcmc = extractchain(mcmc);

%% ROPE
% Compute HDIs:
HDIstruct = structfun(@hdimcmc,mcmc,'UniformOutput',false);
% Define list for recording results:
% Goal: at least one theta greater than ROPE with none below:
goalAchieved.beta4muAboveROPE = HDIstruct.beta4mu(1)>goalgain;
goalAchieved.beta2muAboveROPE = HDIstruct.beta2mu(2)<goalbias;

%% MAE and RMS
[uc,~,subs] = unique([data.s data.c],'rows');
uc = uc(:,2);
mae			= accumarray(subs,abs(data.y-data.x),[],@ mean);
rms			= accumarray(subs,(data.y-data.x).^2,[],@ mean);
rms = sqrt(rms);

% mu			= accumarray(uc+1,mae,[],@mean);
% sd			= accumarray(uc+1+1,mae,[],@std);
h	= ttest(mae(uc==1),mae(~uc));
goalAchieved.MAE = h;
h	= ttest(rms(uc==1),rms(~uc));
goalAchieved.RMS = h;

% keyboard

%%

% Goal: at least one theta greater than ROPE with none below:
b		= mcmc.beta1;
nsub	= size(b,2);
beta1	= NaN(2,nsub);
for subIdx				= 1:nsub
	beta1(:,subIdx)		= hdimcmc(b(:,subIdx));
end
goalAchieved.beta1NarrowHDI = sum((beta1(2,:)-beta1(1,:))<goalbiasnarrow)/nsub;


% Goal: at least one theta greater than ROPE with none below:
b		= mcmc.beta1+mcmc.beta2;
nsub	= size(b,2);
beta	= NaN(2,nsub);
for subIdx				= 1:nsub
	beta(:,subIdx)		= hdimcmc(b(:,subIdx));
end
goalAchieved.biasBestNarrowHDI = sum((beta(2,:)-beta(1,:))<goalbiasnarrow)/nsub;


% Goal: at least one theta greater than ROPE with none below:
b		= mcmc.beta2;
nsub	= size(b,2);
beta2	= NaN(2,nsub);
for subIdx				= 1:nsub
	beta2(:,subIdx)		= hdimcmc(b(:,subIdx));
end
goalAchieved.beta2NarrowHDI = sum((beta2(2,:)-beta2(1,:))<goalbiasnarrow)/nsub;
goalAchieved.beta2AboveROPE = sum(beta2(2,:)<-goalbias)/nsub;

d = data.p.beta2';
d
goalAchieved.beta2hit = sum((beta2(2,:)<-goalbias & d<-goalbias) | (beta2(1,:)>goalbias & d>goalbias))/sum(d>goalbias | d<-goalbias);

% Goal: at least one theta greater than ROPE with none below:
b		= mcmc.beta4;
nsub	= size(b,2);
beta4	= NaN(2,nsub);
for subIdx				= 1:nsub
	beta4(:,subIdx)		= hdimcmc(b(:,subIdx));
end
goalAchieved.beta4NarrowHDI = sum((beta4(2,:)-beta4(1,:))< goalgainnarrow)/nsub;
goalAchieved.beta4AboveROPE = sum(beta4(1,:)>goalgain)/nsub;
d = data.p.beta4';

goalAchieved.beta4hit = sum((beta4(1,:)>goalgain & d>goalgain) | (beta4(2,:)<-goalgain & d<-goalgain))/sum(d>goalgain | d<-goalgain);

% Goal: at least one theta greater than ROPE with none below:
b		= mcmc.beta3;
nsub	= size(b,2);
beta3	= NaN(2,nsub);
for subIdx				= 1:nsub
	beta3(:,subIdx)		= hdimcmc(b(:,subIdx));
end
goalAchieved.beta3NarrowHDI = sum((beta3(2,:)-beta3(1,:))< goalgainnarrow)/nsub;


b		= mcmc.beta3+mcmc.beta4;
nsub	= size(b,2);
beta	= NaN(2,nsub);
for subIdx				= 1:nsub
	beta(:,subIdx)		= hdimcmc(b(:,subIdx));
end
goalAchieved.gainBestNarrowHDI = sum((beta(2,:)-beta(1,:))< goalgainnarrow)/nsub;

%%
% keyboard
% %%
% goalAchieved.beta1AboveROPE = any(beta1(1,:) > nullROPE(2)) & ~any(beta1(2,:) < nullROPE(1) );
% % Goal: all theta's HDI width less than max width:
% goalAchieved.beta1NarrowHDI = all((beta1(2,:)-beta1(1,:))< HDImaxWid );
% % More goals can be inserted here if wanted...


goalAchieved

function [X,C,Y,S] = getdata
clc
cd('/Users/marcw/DATA/Renedata/New');

nsubjects = 14;
Best.Best1 = [];
X = [];
Y = [];
S = [];
C = [];
CIear = [1 1 0 1 0 1 0 0 1 1 1 1 0 1];
CIear = [1 1 0 1 0 1 1 0 1 1 0 1 1 1];
xpassive = linspace(-90,90,33);
xactive = [3 7 11 15 17 19 23 27 31];
xactive = xpassive(xactive);

for subIdx = 1:nsubjects
	strbest		= ['best' num2str(subIdx)];
	strbimodal	= ['bimodal' num2str(subIdx)];
	
	S1	= load(strbest); % upper - median - lower
	S2	= load(strbimodal);
	
	%% Best
	fn	= fieldnames(S1);
	x	= S1.(fn{1});
	
	y = x(:,2);
	x = x(:,1);
	k = dsearchn(xactive',x);
	x = xactive(k)';
	
	[m,n] = size(x);
	x = reshape(x,3,m/3);
	x = mean(x);
	y = reshape(y,3,m/3);
	
	M = y(2,:);
	U = y(1,:)-M;
	L = M-y(3,:);
	
	if ~CIear(subIdx);
		M = -M;
		x = -x;
	end
	sd	= mean([U; L]);
	
	figure(1)
	subplot(3,5,subIdx)
	errorbar(x,M,L,U,'ko','MarkerFaceColor','w')
	hold on
	xlim([-100 100]);
	ylim([-100 100]);
	title(CIear(subIdx))
	
	%% Draw
	nrep = 18;
	nx = numel(x);
	x	= repmat(x,nrep,1);
	M	= repmat(M,nrep,1);
	sd	= repmat(sd,nrep,1);
	
	y = M+sd.*randn(size(M));
	for xIdx = 1:nx
		z = zscore(y(:,xIdx));
		y(:,xIdx) = M(1,xIdx)+z.*sd(1,xIdx);
		
		
	k = dsearchn(xpassive',y(:,xIdx));
	y(:,xIdx) = xpassive(k);	
	
	end
	
	x = x(:);
	y = y(:);
	s = repmat(subIdx,size(x));
	c = ones(size(x));
	
	X = [X x];
	Y = [Y y];
	S = [S s];
	C = [C c];
	
	%% Bimodal
	fn	= fieldnames(S2);
	x	= S2.(fn{1});
	y = x(:,2);
	x = x(:,1);
	k = dsearchn(xactive',x);
	x = xactive(k)';
	
	
	[m,n] = size(x);
	x = reshape(x,3,m/3);
	x = mean(x);
	y = reshape(y,3,m/3);
	M = y(2,:);
	U = y(1,:)-M;
	L = M-y(3,:);
	sd	= mean([U; L]);
	if ~CIear(subIdx);
		M = -M;
		x = -x;
	end
	
	%% Graphics
	subplot(3,5,subIdx)
	errorbar(x,M,L,U,'ko','MarkerFaceColor','k')
	axis square
	box off
	
	xlim([-100 100]);
	ylim([-100 100]);
	unityline;
	title(CIear(subIdx))
	
	%% Draw
	nrep = 18;
	nx = numel(x);
	x	= repmat(x,nrep,1);
	M	= repmat(M,nrep,1);
	sd	= repmat(sd,nrep,1);
	
	y = M+sd.*randn(size(M));
	for xIdx = 1:nx
		z = zscore(y(:,xIdx));
		y(:,xIdx) = M(1,xIdx)+z.*sd(1,xIdx);
		
		
	k = dsearchn(xpassive',y(:,xIdx));
	y(:,xIdx) = xpassive(k);	
		end
	
	x = x(:);
	y = y(:);
	s = repmat(subIdx,size(x));
	c = zeros(size(x));
	
	X = [X x];
	Y = [Y y];
	S = [S s];
	C = [C c];
	
end
X = X(:);
Y = Y(:);
S = S(:);
C = C(:);
% savegraph('gifford','png');
% save('gifford','X','Y','S','C');
%%
[uc,~,subs] = unique([S C],'rows');
uc = uc(:,2);
mae			= accumarray(subs,(Y-X).^2,[],@ mean);
sel = uc==1;
mean(sqrt(mae(sel)))
mean(sqrt(mae(~sel)))
%%
% keyboard


function [xVec,cVec,yVec,sVec,Par] = simdata(groupMCMC,Nsubj,Nrep)

%%  Generate random parameter values for idealized subjects:
b(1).beta		= normrnd(groupMCMC.beta1mu,groupMCMC.beta1sigma,[Nsubj,1]); % bias
b(2).beta		= normrnd(groupMCMC.beta2mu,groupMCMC.beta2sigma,[Nsubj,1]); % bias
b(3).beta		= normrnd(groupMCMC.beta3mu,groupMCMC.beta3sigma,[Nsubj,1]); % bias
b(4).beta		= normrnd(groupMCMC.beta4mu,groupMCMC.beta4sigma,[Nsubj,1]); % bias
Sigma			= groupMCMC.sigma;
Par.beta1 = b(1).beta;
Par.beta2 = b(2).beta;
Par.beta3 = b(3).beta;
Par.beta4 = b(4).beta;
Par.sigma = Sigma;
%% Generate idealized data very close to theta's:
% z					= round(theta*idealNtrl);
% Nrep = 18;
x				= transpose(-75:15:75);
xVec0			= repmat(x,Nrep,1);
xVec1			= xVec0;
xVec			= [xVec0; xVec1];
cVec			= [zeros(size(xVec0));ones(size(xVec1))];
yVec			= NaN(numel(xVec),Nsubj);

%%
for sIdx = 1:Nsubj
	b1 = b(1).beta(sIdx);
	b2 = b(2).beta(sIdx);
	b3 = b(3).beta(sIdx);
	b4 = b(4).beta(sIdx);
	eta = Sigma.*randn(numel(xVec),1);
	yVec(:,sIdx) = b1+b2.*cVec+b3*xVec+b4.*cVec.*xVec + eta;
end

%%
sVec			= repmat(1:Nsubj,[numel(xVec),1]);
xVec			= repmat(xVec,[1,Nsubj]);
cVec			= repmat(cVec,[1,Nsubj]);
yVec			= yVec(:);
sVec			= sVec(:);
xVec			= xVec(:);
cVec			= cVec(:);
