function Rhat = gelmandiag(x,varargin)
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
parNames		= fieldnames(x);
confidence		= keyval('confidence',varargin,0.95);
transform		= keyval('transform',varargin,false);
autoburnin		= keyval('autoburnin',varargin,false);
multivariate	= keyval('multivariate',varargin,true);

nChain			= nchain(x);
if nChain<2
	error('You need at least 2 chains');
end

if autoburnin
	x			= x(round(length(x)/2)+1:end,:);
end
nIter			= niter(x);
nVar			= nvar(x);

if transform
	        x = gelmantransform(x);  % transform the data to improve
	%         normal approximation. Otherwise, confidence limits are incorrectly estimated.
end

%% Mean and variance of each chain
M		= structfun(@(y) mean(y,2),x,'UniformOutput',false); % mean of each chain
S2		= structfun(@(y) var(y,0,2),x,'UniformOutput',false); % variance of each chain

%% Arrays
M		= struct2array(M);
S2		= struct2array(S2);

%% 3. Calculate the within-chain and between-chain variance.
% the mean of the empirical variance within each chain, W
W		= mean(S2);

% B/n is the empirical between-chain variance
B		= nIter*var(M);

%%
if nVar > 1 && multivariate
	%         if (is.R()) {
	%             CW = chol(W)
	
	%             emax = eigen(backsolve(CW, transpose(backsolve(CW, B, transpose = TRUE)),
	%                 transpose = TRUE), symmetric = TRUE, only.values = TRUE).values[1]
	%         }


	%             emax = eigen(qr.solve(W, B), symmetric = FALSE,
	%                 only.values = TRUE).values
	%         mpsrf = sqrt((1 - 1/nIter) + (1 + 1/nVar) * emax/nIter)
	[Q,~]		= qr(W);
	C			= Q'*B;
	[~,S,~]		= svd(C);
	emax		= diag(S)/max(diag(S));
	mpsrf		= sqrt((1 - 1/nIter) + (1 + 1/nVar) * emax/nIter);
else
	mpsrf		= [];
end

%% 4. Calculate the estimated variance of the parameter
% as a weighted sum of the within-chain W and between-chain B variance.
V			= (nIter - 1) * W/nIter + (1 + 1/nChain) * B/nIter;

%% 5. Calculate the potential scale reduction factor Rhat
% see R2.estimate below
% Rhat		= sqrt(V./W); 

%% Distribution
% mean of within-chain mean
muhat		= mean(M);

% variance of within-chain variance
varw		= var(S2);

% variance of between-chain variance
varb		= (2 * B.^2)/(nChain - 1);

% covariance parameters
c1			= cov(S2,M.^2);
c2			= cov(S2,M);
a			= nIter/nChain;
b			= 2*muhat;
covwb		= transpose(a * c1(:) - b*c2(:));

varV		= ( (nIter - 1)^2 * varw + (1 + 1/nChain)^2 * varb +  2 * (nIter - 1) * (1 + 1/nChain) * covwb )/nIter^2;
dfV			= (2 * V.^2)./varV;
dfadj		= (dfV + 3)/(dfV + 1);
Bdf			= nChain - 1;
Wdf			= (2 * W.^2)./varw;
R2.fixed	= (nIter - 1)/nIter;
R2.random	= (1 + 1/nChain) * (1/nIter) * (B./W);
R2.estimate = R2.fixed + R2.random;
upperP		= (1 + confidence)/2;
R2.upper	= R2.fixed + finv(upperP, Bdf, Wdf) .* R2.random;

%% Output
dfadj = 1;
for iVar = 1:nVar
	Rhat.(parNames{iVar}).pointest = sqrt(dfadj * R2.estimate);
	Rhat.(parNames{iVar}).upperCI = sqrt(dfadj * R2.upper);
end

