function postSummary = summarizepost(paramSampleVec,varargin)
% POSTSUMMARY = SUMMARIZEPOST(P)
%
% Obtain summary statistics in structure POSTSUMMARY of posterior P.
%
% Summary statistics include:
% - mean
% - median
% - mode
% - density (estimated via ksdensity)
% - highest-density low and high interval
% - highest-density mass
% - comparison to value
% - ROPE
%
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij

%% input
credMass	= keyval('credMass', varargin); if isempty(credMass), credMass = 0.95; end;
compVal		= keyval('compVal', varargin);
ROPE		= keyval('ROPE', varargin);

% Determine interesting parameters
postSummary.mean	= mean(paramSampleVec);
postSummary.median	= median(paramSampleVec);
[mcmcDensity.y,mcmcDensity.x] = ksdensity(paramSampleVec);
[~,indx]			= max(mcmcDensity.y);
postSummary.mode	= mcmcDensity.x(indx);
HDI					= hdimcmc(paramSampleVec,credMass);
postSummary.hdiMass = credMass;
postSummary.hdiLow	= HDI(1);
postSummary.hdiHigh	= HDI(2);
postSummary.density = mcmcDensity;

%% the comparison value
if ~isempty(compVal)
	postSummary.compVal		= compVal;
	postSummary.pcGTcompVal = sum(paramSampleVec>compVal)/length(paramSampleVec);
end
%% the ROPE
if ~isempty(ROPE)
	pcInROPE				= sum(paramSampleVec>ROPE(1) & paramSampleVec<ROPE(2))/length(paramSampleVec);
	postSummary.ROPElow		= ROPE(1);
	postSummary.ROPEhigh	= ROPE(2);
	postSummary.pcInROPE	= pcInROPE;
end