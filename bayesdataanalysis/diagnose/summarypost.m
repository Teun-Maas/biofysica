function [postSummary,histInfo] = summarypost(paramSampleVec, varargin)
% PLOTPOST(P)
%
% plot MCMC parameter sample vector P
%
% POSTSUMMARY = PLOTPOST(P)
%
% credMass	= keyval('credMass', varargin,0.95);
% showMode	= keyval('showMode',false);
% showCurve	= keyval('showCurve', varargin,false); 
% compVal		= keyval('compVal', varargin);
% ROPE		= keyval('ROPE', varargin);
% yaxt		= keyval('yaxt', varargin);
% ylab		= keyval('ylab', varargin);
% xlab		= keyval('xlab', varargin,'Parameter'); 
% xl			= keyval('xlim', varargin,minmax([compVal; paramSampleVec]')); 
% main		= keyval('main', varargin);
% col			= keyval('col', varargin,[.6 .6 1]);
% breaks		= keyval('breaks', varargin);
%
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij


%% Check arguments
% At least this is better in R
% Override defaults of hist function, if not specified by user
if nargin<1
	paramSampleVec = 1;
end
credMass	= keyval('credMass', varargin,0.95);
compVal		= keyval('compVal', varargin);
ROPE		= keyval('ROPE', varargin);
main		= keyval('main', varargin);
breaks		= keyval('breaks', varargin);
% HDItextPlace = keyval('HDItextPlace', varargin); if isempty(HDItextPlace), HDItextPlace = 0.7; end;
% border		= keyval('border', varargin); if isempty(border), border='w'; end;

%% Determine interesting parameters
postSummary.mean	= mean(paramSampleVec);
postSummary.median	= median(paramSampleVec);
[mcmcDensity.y,mcmcDensity.x] = ksdensity(paramSampleVec);
[~,indx]			= max(mcmcDensity.y);
postSummary.density.x	= mcmcDensity.x;
postSummary.density.y	= mcmcDensity.y;

postSummary.mode	= mcmcDensity.x(indx);
HDI					= hdimcmc(paramSampleVec,credMass);
postSummary.hdiMass = credMass;
postSummary.hdiLow	= HDI(1);
postSummary.hdiHigh	= HDI(2);



if isempty(breaks)
	by=(HDI(2)-HDI(1))/18;
	breaks = unique([min(paramSampleVec):by:max(paramSampleVec) max(paramSampleVec)]);
end
N					= histc(paramSampleVec,breaks);
db					= mean(diff(breaks));
histInfo.N			= N;
histInfo.density	= N./db/sum(N);

%% Display the comparison value.
if ~isempty(compVal)
	cvCol					= [0 .3 0];
	pcgtCompVal				= round(100*sum(paramSampleVec>compVal)/length(paramSampleVec)); % percentage greater than
	pcltCompVal				= 100 - pcgtCompVal; % percentage lower than
	plot([compVal compVal],[0.96*cvHt 0],'k-','Color',cvCol,'LineWidth',2);
	str						= [num2str(pcltCompVal) '% < ' num2str(compVal,3) ' < ' num2str(pcgtCompVal) '%'];
	text(compVal,cvHt,str,'HorizontalAlignment','center','Color',cvCol);
	postSummary.compVal		= compVal;
	postSummary.pcGTcompVal = sum(paramSampleVec>compVal)/length(paramSampleVec);
end
%% Display the ROPE.
if ~isempty(ROPE)
	pcInROPE				= sum(paramSampleVec>ROPE(1) & paramSampleVec<ROPE(2))/length(paramSampleVec);
	postSummary.ROPElow		= ROPE(1);
	postSummary.ROPEhigh	= ROPE(2);
	postSummary.pcInROPE	= pcInROPE;
end


