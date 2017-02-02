function [X,Y,Z] = pa_yaw(X,Y,Z,Angle)
% [X,Y,Z] = PA_YAW(X,Y,Z,YAW);
%
% Rotate YAW (deg) about the Z-axis
%
% See also PA_ROLL, PA_PITCH

% 2007 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
% Compensate for Coordinate Definition
[mx,nx] = size(X);
[my,ny] = size(Y);
[mz,nz] = size(Z);

% Compensate for Coordinate Definition in Matlab
x = X(:);
y = -Z(:);
z = Y(:);

Angle = pa_deg2rad(Angle); % define in radians

%% Define Roll Rotation Matrix for 3D
R = [cos(Angle)  0 sin(Angle)  0;...
     0           1 0           0;...
     -sin(Angle) 0 cos(Angle)  0;...
     0           0 0           1];

%% Define Data Matrix
M = [x y z ones(size(x))]';

%% Matrix Multiply
M = R*M;

%% Revert
M = M';
X = M(:,1);
Y = M(:,3);
Z = -M(:,2);

% Reshape
X = reshape(X,mx,nx);
Y = reshape(Y,my,ny);
Z = reshape(Z,mz,nz);
