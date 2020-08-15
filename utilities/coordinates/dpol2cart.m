function varargout = dpol2cart(A,E,R)
% [X,Y,Z] = DPOL2CART(A,E,R)
%
% transforms corresponding elements of data stored in double-polar (azimuth
% angle A, elevation angle E, and radius R) to Cartesian [X,Y,Z]
% coordinates. The arrays A, E and R must be the same size. Angles A and
% E are in degrees.  
%
% If A is a matrix with three columns (and E and R are not given), they are
% taken to be [A E R]. If A is a matrix with two columns (and E and R are not given), they are
% taken to be [A E] and R is set to unity. 
%
% See also DPOL2CART, DPOL2SPH, SPH2DPOL, DPOL2IPOL, IPOL2DPOL

% AUTHOR: Marc M. van Wanrooij

%% Initialization
if nargin==1 && size(A,2)==2
  E			= A(:,2);
  A			= A(:,1);
end
if nargin<3
    R		= ones(size(A));
end

%% Express angles in radians 
D2R			= pi/180; % radians to degrees
A			= A*D2R;
E			= E*D2R;

%% Convert to Cartesian Coordinates
X			= R*sin(A); % Left-right
Y			= R*sin(E); % Up-Down

signZ		= sign(cos(A).*cos(E)); % front-back
absZ		= abs(sqrt(R.^2-X.^2-Y.^2));
Z			= signZ .* absZ;


%% Output
if nargout == 1 || nargout == 0
    varargout(1)    = {[X Y Z]};
elseif nargout == 3
    varargout(1)    = {X};
    varargout(2)    = {Y};
    varargout(3)    = {Z};
else
    disp('Wrong number of output arguments');
end