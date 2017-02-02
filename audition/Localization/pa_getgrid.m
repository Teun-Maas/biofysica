function [azimuth,elevation,dAz,dEl]=pa_getgrid(nrgridpoints,lngsquare)
% GETGRID Generate grid-positions
%
% GETGRID(NRGRID,RANGE)
% Generate a grid of equidistant NRGRID*NRGID locations, with NRGID
% azimuths and NRGRID elevations. The total range of azimuths or elevations
% will be chosen from RANGE-2*DAZ/DEL deg (default: sqrt(2*90^2) -> DAZ and
% DEL will be -25.4558 deg).
%
% GETGRID is commonly used for GENEXPERIMENT, to create a square grid that
% is subsequently rotated by 45 deg to fit the double-polar coordinate
% system.
%
% Example
% To create a grid of 5 azimuth locations and 5 elevation
% locations, rotate the matrix to fit the double polar system, and plot
% these locations:
%   [azimuth, elevation] = getgrid(5);
%   [azimuth, elevation]= rotategrid(azimuth,elevation,45);
%   cla;
%   plotfart;
%   hold on;
%   h = plot(azimuth,elevation,'ko');
%   set(h,'MarkerFaceColor','w');
%   xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
%
% See also GENEXPERIMENT

% Author MW
% Copyright 2007

if nargin<2
    lngsquare       = sqrt(2*90^2);
end
dAz                 = lngsquare/nrgridpoints;   % deg
dEl                 = lngsquare/nrgridpoints;   % deg
az                  = -(lngsquare/2-dAz/2):dAz:(lngsquare/2-dAz/2);
el                  = -(lngsquare/2-dEl/2):dEl:(lngsquare/2-dEl/2);
[azimuth,elevation] = meshgrid(az,el);