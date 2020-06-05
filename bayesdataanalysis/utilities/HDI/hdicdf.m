function HDI = hdicdf(ICDFname,credMass,tol,varargin)
% HDI = HDIofCDF(ICDFNAME,CREDMASS,TOL)
%
% Arguments:
%   ICDFname is Matlab's name for the inverse cumulative density function
%     of the distribution.
%   credMass is the desired mass of the HDI region.
%   tol is passed to Matlab's optimize function.
% Return value:
%   Highest density iterval (HDI) limits in a vector.
% Example of use: For determining HDI of a beta(30,12) distribution, type
% >> shape1 = 30;
% >> shape2 = 12;
% >> HDIofICDF('beta',shape1,shape2);
%
% ONLY WORKS FOR 2-PARAMETER DISTRIBUTIONS
%
% Original Adapted and corrected from Greg Snow's TeachingDemos package.
% Original in R:	Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
%					A Tutorial with R and BUGS. Academic Press / Elsevier.
% Modified to Matlab code: Marc M. van Wanrooij



%% Check
if nargin<1
	ICDFname = 'beta';
end
if nargin<2
	credMass = 0.95;
end
if nargin<3
	tol = 1e-08;
end
options			= optimset('TolFun',tol);
incredMass		=  1.0 - credMass;

%% Optimize
HDIlowTailPr	= fminbnd(@(lowTailPr) intervalWidth(lowTailPr,ICDFname,credMass,varargin{1},varargin{2}),0,incredMass,options);
HDI				= [icdf(ICDFname,HDIlowTailPr,varargin{1},varargin{2}) icdf(ICDFname,credMass + HDIlowTailPr,varargin{1},varargin{2})];


function iw = intervalWidth(lowTailPr,ICDFname,credMass,varargin)
iw		= icdf(ICDFname,credMass+lowTailPr,varargin{1},varargin{2}) - icdf(ICDFname,lowTailPr,varargin{1},varargin{2});