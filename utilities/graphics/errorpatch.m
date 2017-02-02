function [hpatch,hline] = errorpatch(X,Y,E,col)
% ERRORPATCH(X,Y,E)
%
% plots the graph of vector X vs. vector Y with error patch specified by
% the vector E.
%
% ERRORPATCH(...,'ColorSpec') uses the color specified by the string
% 'ColorSpec'. The color is applied to the data line and error patch, with
% the error patch having an alpha value of 0.4.
%
% [HPATCH,HLINE] = PA_ERRORPATCH(...) returns a vector of patchseries and
% lineseries handles in HPATCH and HLINE, respectively.

% (c) 2011 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

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
if size(E,1)>2
	E   = E(:)';
	if size(E,1)>2
		error('E should be a row vector or 2-row matrix');
	end
end
if length(Y)~=length(X)
	error('Y and X should be the same size');
end
if size(E,2)~=size(X,2)
	error('E and X should be the same size');
end
if nargin<4
	col = 'k';
end

%% remove nans
if size(E,1)>1
	sel		= isnan(X) | isnan(Y) | isnan(E(1,:)) | isnan(E(2,:));
	E		= E(:,~sel);
else
	sel		= isnan(X) | isnan(Y) | isnan(E);
	E		= E(~sel);
end
X		= X(~sel);
Y		= Y(~sel);

%% Create patch
x           = [X fliplr(X)];
if size(E,1)>1
	y           = [E(1,:) fliplr(E(2,:))];
else
	y           = [Y+E fliplr(Y-E)];
end
%% Graph
hpatch           = patch(x,y,col);
alpha(hpatch,0.4);
set(hpatch,'EdgeColor','none');
hold on;
hline = plot(X,Y,'k-'); set(hline,'LineWidth',2,'Color',col);
box on;