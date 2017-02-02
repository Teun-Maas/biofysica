% function azel = rphiazel(R,Phi)
%
%   Coordinate transformation (r,phi) -> (az,el)
%
%   Jeroen Goossens

function [XY] = pa_rphi2azel(R,Phi)

if nargin==1 && size(R,2)==2
  Phi = R(:,2);
  R   = R(:,1);
end

xyz = pa_rphi2xyz(R,Phi);
XY = pa_xyz2azel(xyz(:,1),xyz(:,2),xyz(:,3));
% if nargout>1
% 	X = XY(:,1);
% 	Y = XY(:,2);
% else
% 	
% end

