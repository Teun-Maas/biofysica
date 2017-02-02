function [hpatch,hline] = pa_stairspatch(X,Y,varargin)
% PA_STAIRSPATCH(X,Y)
%
% plots the graph of vector X vs. vector Y with error patch specified by 
% the vector E.
%
% PA_ERRORPATCH(...,'ColorSpec') uses the color specified by the string
% 'ColorSpec'. The color is applied to the data line and error patch, with
% the error patch having an alpha value of 0.3.
%
% [HPATCH,HLINE] = PA_ERRORPATCH(...) returns a vector of patchseries and 
% lineseries handles in HPATCH and HLINE, respectively.

% 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@gmail.com

%% Initialization
% Check whether 
if size(X,1)>1
	X=X(:)';
	if size(X,1)>1
		error('X should be a row vector');
	end
end
if size(Y,1)>1
	Y   = Y(:)';
	if size(Y,1)>1
		error('Y should be a row vector');
	end
end
if length(Y)~=length(X)
	error('Y and X should be the same size');
end

clr        = pa_keyval('Color',varargin);
if isempty(clr)
	clr			= 'k';
end
alfa        = pa_keyval('alpha',varargin);
if isempty(alfa)
	alfa		= .4;
end
%% remove nans
sel		= isnan(X) | isnan(Y);
X		= X(~sel);
Y		= Y(~sel);

%% Create stairs
dx		= mean(diff(X));
dx		= dx-dx/10000;
X		= repmat(X,2,1);
X(1,:)	= X(1,:)-dx/2;
X(2,:)	= X(2,:)+dx/2;
Y		= repmat(Y,2,1);
X		= X(:);
Y		= Y(:);
XY		= [X Y];
XY		= sortrows(XY,1);
X		= XY(:,1)';
Y		= XY(:,2)';

%% Create patch
x       = [X fliplr(X)];
y       = [Y repmat(min(Y),size(Y))];

%% Graph
hpatch           = patch(x,y,clr);
alpha(hpatch,alfa);
set(hpatch,'EdgeColor',clr);
hold off;
box on;