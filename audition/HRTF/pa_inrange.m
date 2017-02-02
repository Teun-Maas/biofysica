function sel = pa_inrange (x, varargin)
% SEL = PA_INRANGE(X,XMIN,XMAX)
%
% Return selection vector SEL indicating which elements of X lie within
% [XMIN,XMAX].
%
% >> pa_inrange([3 5],[4 1; 5 5])
% returns [0 1]

% 2013 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

if nargin == 2
	xmnmx	= varargin{1};
	xmin	= xmnmx(1);
	xmax	= xmnmx(2);
elseif nargin==3
	xmin	= varargin{1};
	xmax	= varargin{2};
end;

sel = x>=xmin & x<=xmax;
