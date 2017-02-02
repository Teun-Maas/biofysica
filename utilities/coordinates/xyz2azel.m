% azel = xyz2azel(X,Y,Z)
%
%  coordinate trazelformation (x,y,z) -> (az, el). If x is a matrix with three
%  columns (and y and z are not given), they are taken to be [x y z].
%

function azel=xyz2azel(x,y,z)

RTD = (180.0 / pi);

if (nargin==1 && size(x,2)==3)
  y=x(:,2);
  z=x(:,3);
  x=x(:,1);
end

azel = zeros(length(x),2);
azel(:,1) = RTD * atan2 (y, sqrt (x.^2 + z.^2));
azel(:,2) = RTD * atan2 (x, sqrt (y.^2 + z.^2));
