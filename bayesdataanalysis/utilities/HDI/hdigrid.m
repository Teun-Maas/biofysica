function HDIinfo = HDIofGrid(probMassVec,credMass)
% [IDX,HDIMASS,HDIHEIGHT] = HDIofGRID(PROBMASSVEC,CREDMASS)
%
% Arguments:
%   probMassVec is a vector of probability masses at each grid point.
%   credMass is the desired mass of the HDI region.
% Return value:
%   A structure with components:
%   indices is a vector of indices that are in the HDI
%   mass is the total mass of the included indices
%   height is the smallest component probability mass in the HDI
% Example of use: For determining HDI of a beta(30,12) distribution
%   approximated on a grid:
%   > probDensityVec = betapdf( linspace(0,1,201),30,12);
%   > probMassVec = probDensityVec/sum(probDensityVec);
%   > HDIinfo = HDIofGrid(probMassVec);
%   > HDIinfo
%
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij


%% Check input
if nargin<2
	credMass = 0.95;
end
sortedProbMass	= sort(probMassVec,'descend');
HDIheightIdx	= find(cumsum(sortedProbMass)>=credMass,1);
HDIheight		= sortedProbMass(HDIheightIdx);
HDImass			= sum(probMassVec(probMassVec>=HDIheight));
idx				= find(probMassVec>=HDIheight);
HDIinfo.idx		= idx;
HDIinfo.height	= HDIheight;
HDIinfo.mass	= HDImass;

