function [AZ,EL] = pa_fart2azel(X,Y,SKY)
% FART2AZEL Transfrom FART to double polar coordinates.
% [AZ,EL] = FART2AZEL(X,Y) transforms corresponding elements of data
%     stored in FART (Fast Auditory Rotating Target) coordinates X 
%    (Theta deg),Y (speakernumber) to double polar coordinates (azimuth AZ and 
%     elevation EL).  The arrays X and Y must be the same size (or either 
%     can be scalar). AZ and EL are returned in deg.
%
%   WHY COORDINATION TRANSFORMATION:
%       FART coordinates *does not* directly reflect Azimuth-Elevation
%       double-polar coordinate system. The elevation of targets is
%       given in speaker number, which should be transformed to deg
%       (each speaker always has the same elevation according to
%       the double polar system irrespective of FART rotation angle
%       Theta).
%       Azimuth AZ can be deduced by applying the formula:
%
%           AZ = arcsin( sin(Theta) * cos(EL) )
%
%       where Theta = FART rotation angle (deg), and EL = elevation angle
%       (deg).
%
%       Note that the double-polar coordinate system by itself does not
%       allow for any front-back distinction. Therefore, elevation is
%       defined as being absolutely larger than 90 deg when coordinates are
%       positioned in the rear hemifield.
%
% See also AZEL2CART, AZEL2POL, CART2AZEL
%
% MarcW 2007

%% Initialization
% Convert to vectors
[M,N]               = size(X);
X                   = X(:);
Y                   = Y(:);

%% Check for old bug (before 19 March 2007)
% Anything below -200 should be recalculated
% sel             = X<-200;
% X(sel)          = -X(sel)-300;

%% Define Front&Rear Speakers, Front&RearHemifields
selfront            = Y<100;                        % Front Speakers Numbers: 1-29
selrear             = Y>100;                        % Rear Speakers Numbers 101-129
selled30            = Y == 30;                      % ET-LED1 Number: 30
selled31            = Y == 31;                      % ET-LED2 Number: 31
selleds             = selled30 | selled31;          % Both ET-LED Numbers: 30&31
selrearhemi1        = (X>90 | X<-90) & selfront;    % Rear hemifield for front speakers
selrearhemi2        = (X<90 & X>-90) & selrear;   % Rear hemifield for rear speakers

%% Speaker Phi 2 Elevation
EL                  = Y;
EL(selfront)        = -55+(Y(selfront)-1)*5;
EL(selrear)         = -57.5+(Y(selrear)-100-1)*5;
EL(selleds)         = -2.5;

%% Theta 2 Azimuth

% Obtain theta coordinates for the E.T.-L.E.D.S.
X(selled30)         = X(selled30)-14;
X(selled31)         = X(selled31)+14;

% And convert
AZ                  = asind(sind(X).*cosd(EL));
% Correct for rear speakers
AZ(selrear)         = -AZ(selrear);

%% Define rear quadrant
% as having absolute Elevations larger than 90 degrees,
% i.e. the front speakers (1-31) will be positioned in the rear when FART has
% moved farther than 90 deg.
% Also, the rear speakers (101-129) will then be positioned in the front
EL(selrearhemi1)    = sign(EL(selrearhemi1)).*(180-abs(EL(selrearhemi1)));
sel                 = selrearhemi1 & EL == 0;
EL(sel)             = -180;
EL(selrearhemi2)    = sign(EL(selrearhemi2)).*(180-abs(EL(selrearhemi2)));

%% Convert back to MTX
AZ                  = reshape(AZ,M,N);
EL                  = reshape(EL,M,N);
