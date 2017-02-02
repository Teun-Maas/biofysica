function varargout = pa_azel2cart(AZ,EL,R)
% AZEL2CART Transform double polar to Cartesian coordinates.
% [X,Y,Z] = AZEL2CART(AZ,EL) transforms corresponding elements of data
%     stored in double polar coordinates azimuth AZ and elevation EL 
%     to Cartesian coordinates X , Y and Z.  The arrays AZ and EL 
%     must be the same size (or all can be scalar, in deg). 
%
%
% See also CART2AZEL, AZEL2POL, FART2AZEL
%
% MarcW 2007
%% Initialization
if nargin==1 && size(Az,2)==2
  EL                = AZ(:,2);
  AZ                = AZ(:,1);
end
AZ                  = pa_deg2rad(AZ);
EL                  = pa_deg2rad(EL);
if nargin<3
    R               = 1;
end

%% Cartesian Coordinates
X                   = R*sin(AZ); % Left-right
Y                   = R*sin(EL); % Up-Down

signZ               = sign(cos(AZ).*cos(EL));
absZ                = abs(sqrt(R.^2-X.^2-Y.^2));
Z                   = signZ .* absZ;


%% Check whether on unit sphere
% margin              = 0.1;
% OnUnitSphere        = (X.^2 + Y.^2 + Z.^2)<(1.0 + margin);
% X(~OnUnitSphere)    = NaN;
% Y(~OnUnitSphere)    = NaN;
% Z(~OnUnitSphere)    = NaN;
% if any(OnUnitSphere)
%     str             = {'AZEL2CART: Beware!';'Azimuth and elevation coordinates are not on Unit Sphere'};
%     disp(char(str));
% end


%% Output
if nargout == 1
    varargout(1)    = {[X Y Z]};
elseif nargout == 3
    varargout(1)    = {X};
    varargout(2)    = {Y};
    varargout(3)    = {Z};
else
    disp('Wrong number of output arguments');
end