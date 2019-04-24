function pearson = pearson_jags(x,y)
%% Correlation Coefficient
chain_globals;

x		= [x y];

% Constants
[n,~]	= size(x);

%% Sampling
% MCMC Parameters
nburnin		= 500; % How Many Burn-in Samples?
nsamples	= 1000;  %How Many Recorded Samples?
nthin		= 1; % How Often is a Sample Recorded?

% Assign Matlab Variables to the Observed Nodes
datastruct = struct('x',x,'n',n);

%% INTIALIZE THE CHAINS.
nChains		= nChainsDefault; % Number of chains to run.
initsStruct = struct([]);
for ii = 1:nChains
	initsStruct(ii).r		= 0;    % because data are standardized
	initsStruct(ii).mu		= zeros(1,2);        % because data are standardized
	initsStruct(ii).lambda	= ones(1,2);  % because data are standardized
end

modelname = which('pearson.txt');
% Use JAGS to Sample
samples = matjags( ...
	datastruct, ...
	modelname, ...
	initsStruct, ...
	'doparallel' , 1, ...
	'nchains', nChains,...
	'nburnin', nburnin,...
	'nsamples', nsamples, ...
	'thin', nthin, ...
	'monitorparams', {'r','mu','sigma'}, ...
	'savejagsoutput',0, ...
	'verbosity' ,0, ...
	'dic',0,...
	'cleanup',1, ...
	'workingdir' , 'tmpjags' );

r = [samples.r];
r = r(:);

mu		= [samples.mu];
[m,n,k] = size(mu);
mu		= reshape(mu,m*n,k);

sigma		= [samples.sigma];
[m,n,k]		= size(sigma);
sigma		= reshape(sigma,m*n,k);

pearson		= struct('r',r','mu',mu,'sigma',sigma);
