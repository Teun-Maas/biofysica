function B = unzscore(varargin)
% B = UNZSCORE(ZETA0,ZETA,MUX,SDX,MUY,SDY)
%
% Undo Z-transformation, as done by ZSCORE, or ZSCORE
%
% To un-zscore variables:
% B = unzscore(z,m,s);
%	or
% B = unzscore(zx,mx,sx,zy,my,sy);
%
% To undo z-transform for linear regresssion parameters:
% B = unzscore(zeta0,zeta,mx,sx); % if only predictor x was z-transformed
%	or
% B = unzscore(zeta0,zeta,mx,sx,my,sy);
%
% B = a structure containing the fields depending on the input:
%	- beta0
%	- beta
%	- x
%	- y

% 2013 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

narg = nargin;
if narg==6 && numel(varargin{1})==1
	narg = narg+1;
end
switch narg
	case 3 	%% Variable
		z = varargin{1};
		m = varargin{2};
		s = varargin{3};
		x = z*s+m;
	case 6
		zx	= varargin{1};
		mx	= varargin{2};
		sx	= varargin{3};
		zy	= varargin{4};
		my	= varargin{5};
		sy	= varargin{6};
		x	= zx.*sx+mx; % can have multiple columns/predictors/covariates
		y	= zy*sy+my; % one dependent variable
	case 4 		%% Parameters for linear regression
		zeta0	= varargin{1};
		zeta	= varargin{2};
		mx		= varargin{3};
		sx		= varargin{4};
		beta0	= zeta0-zeta.*mx./sx;
		beta	= zeta./sx;
	case 7
		%% Parameters for linear regression
		zeta0	= varargin{1};
		zeta	= varargin{2};
		mx		= varargin{3};
		sx		= varargin{4};
		my		= varargin{5};
		sy		= varargin{6};
		beta0	= zeta0.*sy+my-zeta.*sy.*mx./sx;
		beta	= zeta.*sy./sx;
end
%% Output
if exist(beta,'var')
	B.beta0 = beta0;
	B.beta	= beta;
end
if exist(x,'var')
	B.x		= x;
end
if exist(y,'var')
	B.y		= y;
end