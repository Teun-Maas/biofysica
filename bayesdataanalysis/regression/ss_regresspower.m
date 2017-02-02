function ss_regresspower
% % Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa-Power.R
close all; % This closes all of Matlab's graphics windows.
seed = 47405; rng(seed); % Optional, merely for replicability.

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
[xVec,cVec,yVec,sVec,iVec] = getdata;


%% Analyse per subject
us = unique(sVec);
ns = numel(us);
uc = unique(cVec);
nc = numel(uc);
ns = 1;
nc = 1;
for ii = 1:ns
	for jj = 1:nc
		sel = sVec==us(ii) & cVec==uc(jj);
		if sum(sel)
		x = xVec(sel);
		y = yVec(sel);
		figure(ii)
		subplot(1,3,jj)
		plotloc(x,y);
			drawnow

		end
	end
end
return
%% Run Bayesian analysis on idealized data:
mcmc			= hxregjags(xVec,cVec,yVec,sVec);

%%
figure(800)
clf
m = mcmc.beta1mu;
s = mcmc.beta1sigma;
y = m + s.*randn(size(m));
subplot(311)
plotpost(y,'xlim',[-90 90])

subplot(312)
plotpost(m,'xlim',[-90 90])

subplot(313)
plotpost(mcmc.beta1(:),'xlim',[-90 90])

[m,n] = size(mcmc.beta1);
figure(801)
clf
for ii = 1:n
subplot(4,4,ii)
m = mcmc.beta1(:,ii);
% m = m(:,selCI);
m = m(:);
plotpost(m,'xlim',[-80 80])

verline
end

%%

figure(666)
clf
subplot(131)
plotpost(mcmc.beta1mu)

subplot(132)
plotpost(mcmc.beta2mu)

subplot(133)
plot(mcmc.beta1mu)
%%

keyboard
return

%%
close all; 

subplot(221)
plotpost(mcmc.beta1mu)
subplot(222)
plotpost(mcmc.beta2mu)

subplot(223)
plotpost(mcmc.beta3mu)
subplot(224)
plotpost(mcmc.beta4mu)

%%
keyboard


%% Simulation
% Specify sample size for each simulated data set:
Nsubj				= 10; Ntrl = 15;  % 150 trials total
% Specify the number of simulated experiments:
n					= numel(mcmc.beta1(:));
nSimulatedDataSets	= min([10,n]); % An arbitrary large number.
% Run the simulated experiments:
simCount			= 0;
if exist('goalTally','var')
	clear goalTally % in case previously run from here down
end
% goalTally = struct([]);
if ~exist('goalTally','var') 		   % if goalTally does not exist, create it
	goalTally = [];
end

for simIdx = 1:nSimulatedDataSets
	simCount	= simCount+1;
	disp(['\n\n==================== Simulation ' num2str(simCount) ' of ' num2str(nSimulatedDataSets) '====================\n\n']);
	% Generate random beta1/0mu and beta1/0sd  for group distribution:
	GroupPar0		=  [mcmc.beta3mu(simIdx) mcmc.beta3sigma(simIdx) mcmc.beta1mu(simIdx) mcmc.beta1sigma(simIdx) ];
	GroupPar1		=  [mcmc.beta4mu(simIdx) mcmc.beta4sigma(simIdx) mcmc.beta4mu(simIdx) mcmc.beta4sigma(simIdx) ]+GroupPar0;
	Sigma			= mcmc.sigma(simIdx);
	
	[xVec,cVec,yVec,sVec] = simdata(GroupPar0,GroupPar1,Sigma,Nsubj,Ntrl);
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
mcmc		= hregjags(data.x,data.y,data.s);
% Convert coda object to matrix for convenience:
% Specify criteria for goals:
nullROPE	= [0.7 0.9];
HDImaxWid	= 0.2;

% Compute HDIs:
HDIstruct = structfun(@hdimcmc,mcmc,'UniformOutput',false);
% Define list for recording results:
% Goal: at least one theta greater than ROPE with none below:
names	= fieldnames(HDIstruct);
names	= names(2:end); % remove general theta
sel		= strncmp('beta1',names,5); % theta per subject
names	= names(sel); 
beta1 = NaN(2,sum(sel));
for subIdx = 1:sum(sel)
	beta1(:,subIdx) = HDIstruct.(names{subIdx});
end
goalAchieved.beta1AboveROPE = any(beta1(1,:) > nullROPE(2)) & ~any(beta1(2,:) < nullROPE(1) );
% Goal: all theta's HDI width less than max width:
goalAchieved.beta1NarrowHDI = all((beta1(2,:)-beta1(1,:))< HDImaxWid );
% More goals can be inserted here if wanted...






function [X,C,Y,S,I] = getdata


cd('/Users/marcw/DATA/LV/Bimodal patients')

% Bimodal
CIear		= [0 0 0 1 1 1 1 1 1 1 1 0 1 1];
patients	= [1:7 9:12 17:19];
npatients	= numel(patients);

%% Load data
X = []; Y = []; S = []; I = []; C = [];
for ii = 1:npatients
	% Bimodal
	cd(num2str(patients(ii)))
	load bimodal_block1
	SupSac	= pa_supersac(Sac,Stim,2,1);
	if patients(ii) == 2
		load bimodal_block1_2
		supsac = pa_supersac(Sac,Stim);
		SupSac = [SupSac; supsac]; %#ok<*AGROW>
	end
	
	sel		= SupSac(:,30) == 1;
	x		= SupSac(sel,23);
	y		= SupSac(sel,8);
	int		= SupSac(sel,29);
	s		= repmat(patients(ii),size(y));
	c = ones(size(y));
	
	if ~CIear(ii)
		x = -x; y = -y;
	end
	
	X = [X;x];
	Y = [Y;y];
	S = [S;s];
	I = [I;int];
	C = [C;c];
	cd ..
	
	% CI
	cd(num2str(patients(ii)))
	if exist('CIonly.mat','file');
		load CIonly
		SupSac	= pa_supersac(Sac,Stim,2,1);
		
		
		sel		= SupSac(:,30) == 1;
		x		= SupSac(sel,23);
		y		= SupSac(sel,8);
		int		= SupSac(sel,29);
		s		= repmat(patients(ii),size(y));
		c = zeros(size(y));
		if ~CIear(ii)
			x = -x; y = -y;
		end
		X = [X;x];
		Y = [Y;y];
		S = [S;s];
		I = [I;int];
		C = [C;c];
	end
	cd ..
	
	
	% HA
	cd(num2str(patients(ii)))
	if exist('HAonly.mat','file');
		load HAonly
		SupSac	= pa_supersac(Sac,Stim,2,1);
		sel		= SupSac(:,30) == 1;
		x		= SupSac(sel,23);
		y		= SupSac(sel,8);
		int		= SupSac(sel,29);
		s		= repmat(patients(ii),size(y));
		c =repmat(2,size(y));
		if ~CIear(ii)
			x = -x; y = -y;
		end
		X = [X;x];
		Y = [Y;y];
		S = [S;s];
		I = [I;int];
		C = [C;c];
	end
	cd ..
end


sel = ~ismember(S,[3 5]);
X = X(sel);
Y = Y(sel);
S = S(sel);
I = I(sel);
C = C(sel);

