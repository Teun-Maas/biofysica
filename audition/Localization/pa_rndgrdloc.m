function [az,el,theta,phi] = pa_rndgrdloc(nrgridpoints,nrloc)
% [AZ,EL] = PA_RNDGRDLOC(NRGRID,NRLOC)
% Obtain Stimulus Positions THETA (FART azimuth [deg]) and PHI (FART
% speaker [#]). 
%
% GETLOC divides the entire double-polar coordinate system in NRGRID
% rhombi, and then determines NRLOC random locations within each rhombus.
% These random locations are chosen such that the average location falls in
% line with the middle of each rhombus. In total GETLOC will obtain
% NRGRID*NRLOC positions.
%
% See also GENEXPERIMENT, GETPOTLOC, GETGRID, ROTATEGRID, RNDLOC, AZEL2FART
%
% Author: MW

if nargin<1
	nrgridpoints = 5;
end
if nargin<2
	nrloc = 3;
end
%
% Get Potental locations
[Xpot,Ypot]                     = getpotloc;
% Get Grid
[azimuth,elevation]             = getgrid(nrgridpoints);
% Rotate Grid
[PostRotAZ,PostRotEL]           = rotategrid(azimuth,elevation);
% Removing some elevation impossibilities in the PostRotation Matrices
sel                             = PostRotEL<-57.5 | PostRotEL>85;
PostRotEL                       = PostRotEL(~sel);
PostRotAZ                       = PostRotAZ(~sel);
% Get Random Potential Locations in Grid
[az,el]                     = rndloc(Xpot,Ypot,PostRotAZ,PostRotEL,nrloc,nrgridpoints);
% And convert to FART coordinates
[theta,phi]                     = pa_azel2fart(az(:),el(:));



function [Xpot,Ypot] = getpotloc
% Generate all possible stimulus locations in the FART-setup
fixdist             = 8;
y                   = -90:0.5:90;
y                   = 2.5*round(y/2.5);
y                   = unique(y);
x                   = -90:1:90;
[Xpot,Ypot]         = meshgrid(x,y);
Xpot                = Xpot(:);
Ypot                = Ypot(:);
% Remove anything not humanly possible
sel                 = (abs(Xpot)+abs(Ypot))>90;
Xpot                = Xpot(~sel);
Ypot                = Ypot(~sel);
% Remove anything not FARTly possible
sel                 = Ypot>85 | Ypot<-57.5;
Xpot                = Xpot(~sel);
Ypot                = Ypot(~sel);
% Remove anything not compatible with fixation LED
theta               = pa_azel2fart(Xpot,Ypot);
sel                 = (theta>-fixdist & theta<fixdist) | theta>180-fixdist | theta<-180+fixdist;
Xpot                = Xpot(~sel);
Ypot                = Ypot(~sel);


function [azimuth,elevation,dAz,dEl]=getgrid(nrgridpoints,lngsquare)
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

function [PostRotAZ,PostRotEL]=rotategrid(azimuth,elevation,Beta)
% ROTATEGRID Rotate matrix by 45 deg (square becomes rhombus)
%
% {ROTAZ,ROTEL] = ROTATEGRID(AZIMUTH,ELEVATION,BETA)
%
% Example
%   % Create a grid
%   [azimuth, elevation] = getgrid(5);
%   % rotate this grid
%   [azimuth, elevation] = rotategrid(azimuth,elevation,45);
%
% See also: GENEXPERIMENT, GETGRID

% Author MW
% Copyright 2007

% Initialization
[m,n]               = size(azimuth);
azimuth             = azimuth(:);
elevation           = elevation(:);
% Rotation Angle
if nargin<3
    Beta            = pa_deg2rad(45);
else
    Beta            = pa_deg2rad(Beta);
end
% Rotation Matrix
RotMat              = [cos(Beta) -sin(Beta) 0;sin(Beta) cos(Beta) 0;0 0 1];
% PreRotation Vector
PreRot              = [azimuth elevation ones(size(elevation))]';
% PostRotation Matrix
PostRot             = RotMat*PreRot;
% Divide in Azimuth and Elevation
PostRotAZ           = PostRot(1,:)';
PostRotEL           = PostRot(2,:)';
% Reshape Mean Grid Locations
PostRotAZ           = reshape(PostRotAZ,m,n);
PostRotEL           = reshape(PostRotEL,m,n);

function [Xrnd,Yrnd]=rndloc(Xpot,Ypot,PostRotAZ,PostRotEL,nrloc,nrgridpoints)
% Generate random locations within the confines of a grid
% Initialization
PostRotAZ                   = PostRotAZ(:);
PostRotEL                   = PostRotEL(:);
Xrnd                        = repmat(NaN*PostRotAZ,1,nrloc);
Yrnd                        = repmat(NaN*PostRotEL,1,nrloc);
% Loop to create NRLOC-1 random locations and obtain the NRLOCth location
% wich averages with the other random locations to the mean grid location
for i                       = 1:length(PostRotAZ)
    % the NRLOCth location may not faill outside the grid rhombus/diamond
    % and the deviation of the average of random locations from the center
    % of each rhombus should not deviate more than a small specified amount
    abssumrnd               = 100;
    absmean                 = 2;
    while abssumrnd>90/nrgridpoints || absmean>0.8
        % first, look for all potential locations in each rhombus
        abssum              = sum([abs(Xpot-PostRotAZ(i)) abs(Ypot-PostRotEL(i))],2);
        sel                 = abssum<(90/nrgridpoints);
        Xgrid               = Xpot(sel);
        Ygrid               = Ypot(sel);
        % from the indexes of these locations
        nrindx              = length(Xgrid);
        % choose NRLOC-1 random indices
        indx                = round(1+(nrindx-1)*rand(nrloc-1,1));
        % which will define the first NRLOC-1 random locations
        Xrnd(i,1:nrloc-1)   = Xgrid(indx);
        Yrnd(i,1:nrloc-1)   = Ygrid(indx);
        % after which we can define the NRLOCth averaging location
        Xrnd(i,nrloc)       = PostRotAZ(i) - sum(Xrnd(i,1:nrloc-1) - repmat(PostRotAZ(i),1,nrloc-1),2);
        Yrnd(i,nrloc)       = PostRotEL(i) - sum(Yrnd(i,1:nrloc-1) - repmat(PostRotEL(i),1,nrloc-1),2);
        % but since the potential locations are all discrete, we must find
        % the closest possible location for the NRLOCth location
        [mdist,mindistindx] = nanmin(hypot(-Xgrid+Xrnd(i,nrloc),-Ygrid+Yrnd(i,nrloc)));
        Xrnd(i,nrloc)       = Xgrid(mindistindx);
        Yrnd(i,nrloc)       = Ygrid(mindistindx);
        % and check whether this location does not fall outside the rhombus
        abssumrnd           = sum([abs(Xrnd(i,nrloc)-PostRotAZ(i)) abs(Yrnd(i,nrloc)-PostRotEL(i))],2);
        % or does not average the other locations to the center of the
        % rhombus
        absmean             = hypot(mean(Xrnd(i,:)-PostRotAZ(i)),mean(Yrnd(i,:)-PostRotEL(i)));
    end
end
