function [PostRotAZ,PostRotEL] = pa_2drotate(azimuth,elevation,Beta)
% [XR,YR] = PA_2DROTATE(X,Y,BETA);

% Rotate matrices X x Y by BETA deg.
%
% See also: PA_ROTATE, PA_ROLL

% (c) 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
[m,n]               = size(azimuth);
azimuth             = azimuth(:);
elevation           = elevation(:);
% Rotation Angle
if nargin<3
    Beta            = pa_deg2rad(45);
else
    Beta            = pa_deg2rad(Beta);
end

%% Rotation Matrix
RotMat              = [cos(Beta) -sin(Beta) 0;sin(Beta) cos(Beta) 0;0 0 1];

%% PreRotation Vector
PreRot              = [azimuth elevation ones(size(elevation))]';

%% PostRotation Matrix
PostRot             = RotMat*PreRot;

%% Divide in Azimuth and Elevation
PostRotAZ           = PostRot(1,:)';
PostRotEL           = PostRot(2,:)';

%% Reshape Mean Grid Locations
PostRotAZ           = reshape(PostRotAZ,m,n);
PostRotEL           = reshape(PostRotEL,m,n);