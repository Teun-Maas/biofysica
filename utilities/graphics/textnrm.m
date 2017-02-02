function [h,Xtxt,Ytxt] = pa_text (x,y,string,varargin)
% PA_TXT(X,Y,STRING,<PROPERTIES>)
%
%    Write a STRING in the current axis, where X and Y are given in
%    normalized axis coordinates. 
%
% See also TEXT, AXIS

% (c) 2011-04-27 Marc van Wanrooij 

%% Get axis values
ax		= axis;
xlim	= ax([1 2]);
ylim	= ax([3 4]);

%% Absolute x- and y-location of text
Xtxt = xlim(1) + x*(xlim(2)-xlim(1));
Ytxt = ylim(1) + y*(ylim(2)-ylim(1));

%% Print to figure
if nargin == 3 % no properties defined
  h = text (Xtxt,Ytxt,string);
else % evaluate properties
  h = text(Xtxt,Ytxt,string,varargin{:});
end;

