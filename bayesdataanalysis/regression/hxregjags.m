function samples = hxregjags(x,c,y,s,varargin)
% MCMC = HXREGJAGS(X,C,Y,S)
%
% Bayesian hierarchical regression with multiplicative interaction:
% Y = beta1 + beta2*X + beta3*X + beta4*C*X
%
% for every subject S, with beta-values drawn from a hierarchical group distribution
% with mean betaNmu and betaNsigma. 
%
% if C denotes condition, coded by 0 and 1, beta2 and beta4 denote
% difference in offset and slope between condition 1 and 0, while beta1 and
% beta3 are offset and slope for condition 0.
%
% MCMC is a structure that contains all beta samples found through MCMC.
% 
% You need to install JAGS and MATJAGS
%	http://mcmc-jags.sourceforge.net/
%	-  http://psiexp.ss.uci.edu/research/programs_data/jags/ and/or https://github.com/msteyvers/matjags
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij
%
% See also REGJAGS, MATJAGS, PLOSTPOST

%% Initialization
chain_globals;
xName			= keyval('xName',varargin,'x');
yName			= keyval('yName',varargin,'y');
sName			= keyval('sName',varargin,'s');
cName			= keyval('cName',varargin,'c');
numSavedSteps	= keyval('numSavedSteps',varargin,3000);
thinSteps		= keyval('thinSteps',varargin,1);
saveName		= keyval('saveName',varargin,'HierMultLinRegressData-Jags-');
nChains			= keyval('nChains',varargin,nChainsDefault);

%% Load data
data			= struct(xName,x,cName,c,yName,y,sName,s);
% fileNameRoot	= 'HierMultLinRegressData-Jags-' ;

%% Generate the MCMC chain:
samples = genMCMC('data',data,'xName',xName,'yName',yName,'sName',sName,'cName',cName,...
	'thinSteps',thinSteps,'saveName',saveName,'numSavedSteps',numSavedSteps,'nChains',nChains);

function samples = genMCMC(varargin)

chain_globals;
data			= keyval('data',varargin);
xName			= keyval('xName',varargin);
yName			= keyval('yName',varargin);
sName			= keyval('sName',varargin);
cName			= keyval('cName',varargin);
numSavedSteps	= keyval('numSavedSteps',varargin);
thinSteps		= keyval('thinSteps',varargin);
saveName		= keyval('saveName',varargin);
nChains			= keyval('nChains',varargin);

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
burnInSteps		= 1000;			% Number of steps to 'burn-in' the samplers.
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
	'doparallel', doparallel, ...
	'nchains', nChains,...
	'nburnin', burnInSteps,...
	'nsamples', nIter, ...
	'thin', thinSteps, ...
	'monitorparams', parameters, ...
	'savejagsoutput' , 1 , ...
	'verbosity', 1 , ...
	'dic',0,...
	'cleanup', 0);

if ~isempty(saveName)
	save([saveName 'Mcmc'],'samples');
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