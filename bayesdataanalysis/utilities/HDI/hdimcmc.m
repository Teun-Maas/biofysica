function HDIlim = hdimcmc(sampleVec,credMass)
% Computes highest density interval from a sample of representative values,
%   estimated as shortest credible interval.
% Arguments:
%   sampleVec
%     is a vector of representative values from a probability distribution.
%   credMass
%     is a scalar between 0 and 1, indicating the mass within the credible
%     interval that is to be estimated.
% Value:
%   HDIlim is a vector containing the limits of the HDI
%
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij


%% Check arguments
if nargin<2
	credMass = 0.95;
end

%%
sortedPts	= sort(sampleVec);
ciIdxInc	= floor(credMass*length(sortedPts)); % number of samples included in CI
nCIs		= length(sortedPts) - ciIdxInc; % number of samples NOT included
ciWidth		= zeros(1,nCIs);
for ii = 1:nCIs
	ciWidth(ii) = sortedPts(ii+ciIdxInc) - sortedPts(ii); % determine credible interval width
end
[~,indx]	= min(ciWidth); % the HDI = shortest credible interval
HDImin		= sortedPts(indx);
HDImax		= sortedPts(indx+ciIdxInc);
HDIlim		= [HDImin , HDImax];
