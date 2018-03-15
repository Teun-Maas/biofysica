function [coefCat,coefMet,coef] = regcat(y,xCat,xMet)
% [COEFNOM,COEFMET,COEF] = REGCAT(Y,XCAT,XMET)
%
% Regression Y on XMET with categorical variables XCAT
%
% If you have the Statistics toolbox, you might prefer:
% >> dat			= dataset(y,xCat,xMet);
% >> dat.xCat		= nominal(dat.xCat);
% >> fit			= LinearModel.fit(dat,'y~xCat+xMet');
% 
% Or if you are using R:
% > xCat	= factor(xCat)
% > lmInfo	= lm(y~xCat+xMet)
%
% See also:
% http://blog.yhathq.com/posts/r-lm-summary.html
% http://www.psychstat.missouristate.edu/multiboo       k/mlt08.htm
% http://www.mathworks.nl/help/stats/group-comparisons-using-categorical-arrays.html

% 2013 Marc van Wanrooij
% e-mail: marcvanwanrooij@gmail.com

%% Check
if nargin<2
% 	xCat	= [0 0 1 1 2 3]';
	xCat = {'Thuis','Thuis','Oost','Oost','West','Best'};
end
if nargin<2
	xMet	= [-1 2 2 3 4 8]';
end
if nargin<1
	y		= [0 1 2 3 4 5]';
end
n				= size(y,1); % number of data points
nMet			= size(xMet,2); % number of metric variables
% nCat			= size(xCat,2); % number of categorical variables

%% dummy variables
[uxCat,~,ic]	= unique(xCat,'stable');
nLevels			= numel(uxCat); % number of levels
dCat			= zeros(n,nLevels);
for ii = 1:n
	dCat(ii,ic(ii)) = 1;
end

%% Combine the dummy and metric variables
x			= [ones(n,1) dCat(:,2:nLevels) xMet];% remove the first dummy variable to avoid overfitting

%% Regress
% Statistics toolbox: 
% x			= [dCat(:,2:nLevels) xMet];% remove the first dummy variable to avoid overfitting
% b			= regstats(y,x,'linear','beta');
coef		= x\y;
coefMet		= coef(nLevels+1:nLevels+nMet);

%% Sum the deflections to 0
% Now, all deflections are described with respect to the first categorical variable
% see also:
% Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
% A Tutorial with R and BUGS. Academic Press / Elsevier.
coefCat		= coef(1:nLevels);
coefCat		= [coefCat(1); coefCat(1)+coefCat(2:end)];
coefCat0	= mean(coefCat);
coefCat		= coefCat-coefCat0;

