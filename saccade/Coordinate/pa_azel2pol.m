function [R,PHI,Z]=pa_azel2pol(AZ,EL)
%AZEL2POL Transform Double Polar to Polar coordinates
% RPHI = AZEL2POL(AZ,EL) transforms corresponding elements of data
%    stored in double polar coordinates AZ,EL to polar coordinates (angle
%    PHI and radius R).  The arrays AZ and EL must be the same size (or
%    either can be scalar). PHI is returned in deg. 
% 
%
% See also AZEL2CART, CART2POL, POL2CART, HYPOT, FART2AZEL
%   MarcW 2007

if nargin==1 && size(AZ,2)==2
  EL        = AZ(:,2);
  AZ        = AZ(:,1);
end
[X,Y,Z]     = pa_azel2cart(AZ,EL);
[PHI,R,Z]   = cart2pol(X,Y,Z);
R           = pa_rad2deg(atan2(R,Z));
PHI         = pa_rad2deg(PHI);
sel         = PHI<0;
PHI(sel)    = PHI(sel)+360;
