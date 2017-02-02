function regresspower
% % Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa-Power.R
close all; % This closes all of Matlab's graphics windows.
seed = 47405; rng(seed); % Optional, merely for replicability.

% Load the functions genMCMC, smryMCMC, and plotMCMC:
% (This also sources DBDA2E-utilities.R)
% source('Jags-Ydich-XnomSsubj-MbinomBetaOmegaKappa.R')

%  Specify idealized hypothesis:
idealGroupPar0		= [0.3 0.1 0 5]; % ideal group parameters, Cond 0 & Cond 1 Gain Mean & SD, Bias Mean & SD
idealGroupPar1		= [0.6 0.1 0 5]; % ideal group parameters, Cond 0 & Cond 1 Gain Mean & SD, Bias Mean & SD
idealNsubj			= 20; % more subjects => higher confidence in hypothesis
idealNtrl			= 20; % more trials => higher confidence in hypothesis
idealSigma			= 5; % does not vary per subject

[xVec,cVec,yVec,sVec] = simdata(idealGroupPar0,idealGroupPar1,idealSigma,idealNsubj,idealNtrl);

%% Run Bayesian analysis on idealized data:
mcmc			= hxregjags(xVec,cVec,yVec,sVec);

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

function samples = hxregjags(x,c,y,s)
% MCMC = HREGJAGS(X,Y,S)
%
% Example for Jags-Ymet-XmetSsubj-MrobustHier.R
%% Load data
xName			= 'X' ; yName = 'Y' ; sName='S'; cName='C';
myData			= struct(xName,x,cName,c,yName,y,sName,s);
fileNameRoot	= 'HierMultLinRegressData-Jags-' ;

%% Generate the MCMC chain:
samples = genMCMC('data',myData,'xName',xName,'yName',yName,'sName',sName,'cName',cName,...
	'numSavedSteps',1000,'thinSteps',1,'saveName',fileNameRoot);

function samples = genMCMC(varargin)

chain_globals;
data			= keyval('data',varargin);
xName			= keyval('xName',varargin,'x');
yName			= keyval('yName',varargin,'y');
sName			= keyval('sName',varargin,'s');
cName			= keyval('cName',varargin,'c');
numSavedSteps	= keyval('numSavedSteps',varargin,1000);
thinSteps		= keyval('thinSteps',varargin,1);
saveName		= keyval('saveName',varargin);
nChains			= keyval('nChains',varargin,nChainsDefault);

%% THE DATA.
y = data.(yName);
x = data.(xName);
c = data.(cName);

% Convert sName to consecutive integers:
[~,~,s] = unique(data.(sName));
ns		= max(s);


%Ntotal = length(y)
% Specify the data in a list, for later shipment to JAGS:
dataStruct = struct('x',x,'c',c,'y',y,'s',s,'cx',c.*x,'Nsubj',ns,'Ntotal',length(y),'Nx',size(x,2));


writemodel;


%% RUN THE CHAINS
parameters = {	'beta1','beta1mu','beta1sigma',...
				'beta2','beta2mu','beta2sigma',...
				'beta3','beta3mu','beta3sigma',...
				'beta4','beta4mu','beta4sigma',...
				'sigma','nu'};
burnInSteps		= 100;			% Number of steps to 'burn-in' the samplers.
nIter			= ceil((numSavedSteps*thinSteps )/nChains); % Steps per chain.
if strcmp(runjagsMethodDefault,'parallel');
	doparallel		= 1; % do use parallelization
else
	doparallel		= 0; % do not use parallelization
end
modelname = fullfile(pwd, 'model.txt');

initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).beta1			= zeros(ns,1);    % because data are standardized
	initsStruct(ii).beta2			= zeros(ns,1);        % because data are standardized
	initsStruct(ii).beta3			= zeros(ns,1);        % because data are standardized
	initsStruct(ii).beta4			= zeros(ns,1);        % because data are standardized
	initsStruct(ii).beta1mu		= 0;        % because data are standardized
	initsStruct(ii).beta2mu		= 0;        % because data are standardized
	initsStruct(ii).beta3mu		= 0;        % because data are standardized
	initsStruct(ii).beta4mu		= 0;        % because data are standardized
	initsStruct(ii).beta1sigma		= 1;        % because data are standardized
	initsStruct(ii).beta2sigma		= 1;        % because data are standardized
	initsStruct(ii).beta3sigma		= 1;        % because data are standardized
	initsStruct(ii).beta4sigma		= 1;        % because data are standardized
	initsStruct(ii).sigma			= 1;        % because data are standardized
	initsStruct(ii).nuMinusOne		= 1;        % because data are standardized
end

% Use JAGS to Sample
samples = matjags( ...
	dataStruct, ...
	modelname, ...
	initsStruct, ...
	'doparallel' , doparallel, ...
	'nchains', nChains,...
	'nburnin', burnInSteps,...
	'nsamples', nIter, ...
	'thin', thinSteps, ...
	'monitorparams', parameters, ...
	'savejagsoutput' , 1 , ...
	'verbosity' , 1 , ...
	'dic',0,...
	'cleanup' , 0);

if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
end

parameterNames = fieldnames(samples); % get all parameter names
for ii = 1:numel(parameterNames)
	n = size(samples.(parameterNames{ii}));
	
	switch numel(n)
		case 2
			samples.(parameterNames{ii}) = samples.(parameterNames{ii})(:);
		case 3
			samples.(parameterNames{ii}) = reshape(samples.(parameterNames{ii}),n(1)*n(2),n(3));
			
	end
end

function writemodel
str = ['# Specify the model:\r\n',...
	'model {\r\n',...
	'\tfor ( i in 1:Ntotal ) {\r\n',...
	'\t\ty[i] ~ dt( beta1[s[i]] + beta2[s[i]] * c[i] + beta3[s[i]] * x[i] + beta4[s[i]] * cx[i] , 1/sigma^2 , nu )\r\n',...
	'\t}\r\n',...
	'\tfor ( j in 1:Nsubj ) {\r\n',...
	'\t\tbeta1[j] ~ dnorm(beta1mu , 1/(beta1sigma)^2 )\r\n',...
	'\t\tbeta2[j] ~ dnorm(beta2mu , 1/(beta2sigma)^2 )\r\n',...
	'\t\tbeta3[j] ~ dnorm(beta3mu , 1/(beta3sigma)^2 )\r\n',...
	'\t\tbeta4[j] ~ dnorm(beta4mu , 1/(beta4sigma)^2 )\r\n',...
	'\t}\r\n',...
	'\t# Priors vague on standardized scale:\r\n',...
	'\tbeta1mu ~ dnorm( 0 , 1/(10)^2 )\r\n',...
	'\tbeta2mu ~ dnorm( 0 , 1/(10)^2 )\r\n',...
	'\tbeta3mu ~ dnorm( 0 , 1/(10)^2 )\r\n',...
	'\tbeta4mu ~ dnorm( 0 , 1/(10)^2 )\r\n',...
	'\tbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'\tbeta2sigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'\tbeta3sigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'\tbeta4sigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'\tsigma ~ dunif( 1.0E-3 , 1.0E+3 )\r\n',...
	'\tnu <- nuMinusOne+1\r\n',...
	'\tnuMinusOne ~ dexp(1/29.0)\r\n',...
	'}\r\n',...
	];

% Write the modelString to a file, using Matlab commands:
fid			= fopen('model.txt','w');
fprintf(fid,str);
fclose(fid);

function [xVec,cVec,yVec,sVec] = simdata(GroupPar0,GroupPar1,Sigma,Nsubj,Ntrl)

%%  Generate random parameter values for idealized subjects:
b(1).beta		= normrnd(GroupPar0(3),GroupPar0(4),[Nsubj,1]); % bias
b(2).beta		= normrnd(GroupPar1(3),GroupPar1(4),[Nsubj,1])-b(1).beta; % delta bias
b(3).beta		= normrnd(GroupPar0(1),GroupPar0(2),[Nsubj,1]); % gain
b(4).beta		= normrnd(GroupPar1(1),GroupPar1(2),[Nsubj,1])-b(3).beta; % gain

%% Generate idealized data very close to theta's:
% z					= round(theta*idealNtrl);
xVec0			= transpose(round(linspace(-90,90,Ntrl)/15)*15); % condition 0 = poor
xVec1			= transpose(round(linspace(-90,90,Ntrl)/15)*15); % condition 1 = good
xVec			= [xVec0; xVec1];
cVec			= [zeros(size(xVec0));ones(size(xVec1))];
yVec			= NaN(Ntrl*2,Nsubj);
for sIdx = 1:Nsubj
	yVec(:,sIdx) = b(1).beta(sIdx)+b(2).beta(sIdx).*cVec+b(3).beta(sIdx)*xVec+b(4).beta(sIdx).*cVec.*xVec + Sigma.*randn(Ntrl*2,1);
end
sVec			= repmat(1:Nsubj,[Ntrl*2,1]);
xVec			= repmat(xVec,[1,Nsubj]);
cVec			= repmat(cVec,[1,Nsubj]);
yVec			= yVec(:);
sVec			= sVec(:);
xVec			= xVec(:);
cVec			= cVec(:);
