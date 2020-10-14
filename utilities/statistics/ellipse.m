function [MU,SD,A,x,y] = ellipse(x,y,varargin)
% [MU,SD,A] = ELLIPSE(X,Y)
%
%  Determine main axes (direction and std) of 2-dimensional matrix [X,Y].
%
%
% ELLIPSE(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
%	'outlier'	- removes outliers, that are 'outlier'*standard deviation
%	away from the cardinal axis means.
%
% See also ELLIPSEPLOT, EIG, COV

% 2011  Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
out         = keyval('outlier',varargin);
if isempty(out)
	out			= []; % (Hz)
end

%% Eigen-values for covariance-matrix
[Veig,Deig]     = eig(cov(x,y)); % Veig gives main axes
A               = bf_rad2deg(atan2(Veig(2),Veig(1))); % angle

%% Delete response >3SD
if ~isempty(out)
	SD			= sqrt(diag(Deig)); % diagonal Deig = variance in the 2 main axes
	xmu			= mean(x);
	ymu			= mean(y);
	
	[Xr,Yr]		= rotate2d(x,y,-A);
	seld		= abs(Xr-mean(Xr))<out*SD(2) & abs(Yr-mean(Yr))<out*SD(1);
	x			= x(seld);
	y			= y(seld);
	
	%% Eigen-values for covariance-matrix
	[Veig,Deig] = eig(cov(x,y)); % Veig gives main axes
	A			= rad2deg(atan2(Veig(2),Veig(1))); % angle
end

%% diagonal Deig = variance in the 2 main axes
SD              = sqrt(diag(Deig));
xmu             = mean(x);
ymu             = mean(y);
MU              = [xmu ymu];

