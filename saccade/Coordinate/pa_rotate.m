function [X,Y,Z] = pa_rotate(X,Y,Z,Roll,Pitch,Yaw)
% [X,Y,Z] = PA_ROTATE(X,Y,Z,ROLL,PITCH,YAW);
%
% Rotate ROLL (deg) about the X-axis, PITCH deg about the Y-axis and YAW
% deg about the Z-axis
%
% See also PA_YAW, PA_ROLL, PA_PITCH

% 2007 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

% Compensate for Coordinate Definition
[mx,nx] = size(X);
[my,ny] = size(Y);
[mz,nz] = size(Z);

% Compensate for Coordinate Definition in Matlab
Xp = X(:);
Yp = -Z(:);
Zp = Y(:);

% Point
P = [Xp Yp Zp]';

% Define Quaternion
Roll    = deg2rad(Roll);
Pitch   = deg2rad(Pitch);
Yaw     = deg2rad(Yaw);

s		= cos(Yaw/2)*cos(Pitch/2)*cos(Roll/2)+sin(Yaw/2)*sin(Pitch/2)*sin(Roll/2);
x		= cos(Yaw/2)*sin(Pitch/2)*cos(Roll/2)+sin(Yaw/2)*cos(Pitch/2)*sin(Roll/2);
y		= sin(Yaw/2)*cos(Pitch/2)*cos(Roll/2)-cos(Yaw/2)*sin(Pitch/2)*sin(Roll/2);
z		= cos(Yaw/2)*cos(Pitch/2)*sin(Roll/2)-sin(Yaw/2)*sin(Pitch/2)*cos(Roll/2);

% Define Matrix
M(1,1)	= 1 -    2*(y^2  +   z^2);
M(1,2)	=        2*(x*y  -   s*z);
M(1,3)	=        2*(x*z  +   s*y);
M(2,1)	=        2*(x*y  +   s*z);
M(2,2)	= 1 -    2*(x^2  +   z^2);
M(2,3)	=        2*(y*z  -   s*x);
M(3,1)	=        2*(x*z  -   s*y);
M(3,2)	=        2*(y*z  +   s*x);
M(3,3)	= 1 -    2*(x^2  +   y^2);
P		= M*P;

%% Revert
P = P';
X = P(:,1);
Y = P(:,3);
Z = -P(:,2);

% Reshape
X = reshape(X,mx,nx);
Y = reshape(Y,my,ny);
Z = reshape(Z,mz,nz);
