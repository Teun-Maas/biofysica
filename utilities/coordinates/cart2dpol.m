function [A,E,R] = cart2dpol(X,Y,Z)
% [A,E,R] = CART2DPOL(X,Y,Z)
%
% transforms corresponding elements of data stored in Cartesian [X,Y,Z] to
% double-polar coordinates  (azimuth angle A, elevation angle E, and radius
% R). The arrays X, Y and Z must be the same size. Angles A and E are
% returned in degrees.
%
% If X is a matrix with three columns (and Y and Z are not given), they are
% taken to be [X Y Z]. 
%
% See also DPOL2CART, DPOL2SPH, SPH2DPOL, DPOL2IPOL, IPOL2DPOL


% AUTHOR: Marc M. van Wanrooij

%% Initialization
if nargin==1
  Y = X(:,2);
  Z = X(:,3);
  X = X(:,1);
end
R2D		= 180/pi; % radians to degrees

%% Coordinate transformation
A		= R2D * atan2(Y, sqrt (X.^2 + Z.^2));
E		= R2D * atan2(X, sqrt (Y.^2 + Z.^2));
sel		= z<0 & y>=0;
E(sel)	= -180-E(sel);
sel		= z<0 & y<0;
E(sel)	= 180-E(sel);
R		= sqrt(sum([X Y Z].^2,2));

%% Output
if nargout==1
	A = [A E R];
end