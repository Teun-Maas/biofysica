% function xyz = rphi2xyz(R,Phi)
%
%   Coordinate transformation (r,phi) -> (x,y,z)
%
%   Jeroen Goossens

function ans=rphi2xyz(R,Phi)

if nargin==1 & size(R,2)==2
  Phi = R(:,2);
  R   = R(:,1);
end

ans = zeros(length(R),3);
DTR = (pi / 180.0); 

ans(:,1) = sin(DTR*R) .* sin(DTR*Phi);
ans(:,2) = sin(DTR*R) .* cos(DTR*Phi);
ans(:,3) = cos(DTR*R);
