function y = pa_asin(X,Par)
% ASINFUN   Function of the arcsine a*asin(b*(X+c))+d.
%   ASINFUN is used by the function FITASIN.
a = Par(1);
b = Par(2);
c = Par(3);
d = Par(4);
y = a*asin(b*(X+c))+d;


