function [Rhat,Rhat95] = shrinkfactor(mcmcStruct,varargin)
% RHAT = SHRINKFACTOR(MCMC,'paramName',PARAMNAME)
%
% See: http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf
%
%
% Steps (for each parameter):
% 1. Run m ? 2 chains of length 2n from overdispersed starting
% values.
% 2. Discard the first n draws in each chain.
% 3. Calculate the within-chain and between-chain variance.
% 4. Calculate the estimated variance of the parameter as a
% weighted sum of the within-chain and between-chain variance.
% 5. Calculate the potential scale reduction factor.
%
% From gelman.diag in R:
%
% There are two ways to estimate the variance of the stationary
% distribution: the mean of the empirical variance within each chain, W,
% and the empirical variance from all chains combined, which can be
% expressed as
%
% sigma.hat^2 = (n-1)W/n + B/n
%
% where n is the number of iterations and B/n is the empirical
% between-chain variance.
%
% If the chains have converged, then both estimates are unbiased. Otherwise
% the first method will underestimate the variance, since the individual
% chains have not had time to range all over the stationary distribution,
% and the second method will overestimate the variance, since the starting
% points were chosen to be overdispersed.
%
% The convergence diagnostic is based on the assumption that the target
% distribution is normal. 
% A Bayesian credible interval can be constructed
% using a t-distribution with mean
%
% mu.hat = Sample mean of all chains combined
%
% and variance
%
% V.hat=sigma.hat2 + B/(mn)
%
% and degrees of freedom estimated by the method of moments
%
% d = 2*V.hat^2/Var(V.hat)
%
% Use of the t-distribution accounts for the fact that the mean and
% variance of the posterior distribution are estimated.
%
% The convergence diagnostic itself is
%
% R=sqrt((d+3) V.hat /((d+1)W)
%
% Values substantially above 1 indicate lack of convergence. If the chains
% have not converged, Bayesian credible intervals based on the
% t-distribution are too wide, and have the potential to shrink by this
% factor if the MCMC run is continued.
%
% RHAT = GELMANDIAG(X,NAME,VALUE)
%
% Gelman and Rubin's ( Gelman-Rubin-Brooks) convergence diagnostic
% also termed:
% - potential scale reduction factor
% - shrink factor
%
% References
%
% Gelman, A., Carlin, J.B., Stern, H.S., Duncon, D.B., Vehtari A. and Rubin
% D.B. (20014) Bayesian Data Analysis. 3E. CRC Press, pages 284-286
%
% Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
%
% Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. Journal of Computational and Graphical Statistics, 7, 434-455.
%
% See also SHRINKFACTOR, PLOTGRBMCMC
%
% Converted from R gelman.diag to Matlab

%% Initialization
parNames	= fieldnames(mcmcStruct);
parName		=  keyval('parName',varargin,parNames{1});
binwidth	= keyval('binwidth',varargin,0);
binsize		= keyval('binsize',varargin,50);

%%
x			= mcmcStruct.(parName);
x			= transpose(x);
n			= length(x);

switch binwidth
	case 0
		[Rhat,Rhat95] = getshrink(x);		
	otherwise
		Rhat		= NaN(n-binsize,1);
		Rhat95		= Rhat;
		for ii = (1:n-binsize)*binwidth
			idx = 1:binsize+ii;
			[Rhat(ii),Rhat95(ii)] = getshrink(x(idx,:));
		end
end

function [Rhat,Rhat95] = getshrink(x)
% Steps (for each parameter):
% 1. Run m ? 2 chains of length 2n from overdispersed starting
% values.
% 2. Discard the first n draws in each chain.
% 3. Calculate the within-chain and between-chain variance.
% 4. Calculate the estimated variance of the parameter as a
% weighted sum of the within-chain and between-chain variance.
% 5. Calculate the potential scale reduction factor.
%


nIter	= niter(x);
nChain	= nchain(x);
%% Mean and variance of each chain
M		= mean(x); % mean of each chain
S2		= var(x); % variance of each chain

%% 3. Calculate the within-chain and between-chain variance.
% the mean of the empirical variance within each chain, W
W		= mean(S2);

% B/n is the empirical between-chain variance
B		= nIter*var(M);

%% 4. Calculate the estimated variance of the parameter
% as a weighted sum of the within-chain W and between-chain B variance.
% V			= (nIter - 1) * W/nIter + (1 + 1/nChain) * B/nIter;

%% 5. Calculate the potential scale reduction factor Rhat
% see R2.estimate below
% Rhat		= sqrt(V./W); 

%% Distribution
% mean of within-chain mean

% variance of within-chain variance
varw		= var(S2);

% variance of between-chain variance

Bdf			= nChain - 1;
Wdf			= (2 * W.^2)./varw;
R2.fixed	= (nIter - 1)/nIter;
R2.random	= (1 + 1/nChain) * (1/nIter) * (B./W);
R2.estimate = R2.fixed + R2.random;
upperP		= (1 + 0.95)/2;
R2.upper	= R2.fixed + finv(upperP, Bdf, Wdf) .* R2.random;

Rhat		= sqrt(R2.estimate);
Rhat95		= sqrt(R2.upper);

