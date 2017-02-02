function [Theta,Phi] = azel2fart(AZ,EL)
% AZEL2FART Transfrom double polar to FART coordinates.
% [Theta,Phi] = AZEl2FART(Theta,Phi) transforms corresponding elements of
% data stored in double polar (azimuth AZ and elevation EL) coordinates to
% the nearest FART coordinates Theta (Theta deg), and Phi (speakernumber).
% The arrays AZ and EL must be the same size (or either can be scalar).
%
%
% Note that the double polar coordinate system is NOT defined for rear
% sound sources. This script "solves" that by assuming the convention that large elevations
% (abs(EL)>90) are from the rear, while azimuth is only defined for the
% range -90 to 90 deg.
%
% See also FART2AZEL, AZEL2CART, AZEL2POL, CART2AZEL

% (c)  21 February 2008
% Author: dr. M.M. van Wanrooij


%% Convert to vectors
[M,N]               = size(AZ);
AZ                  = AZ(:);
EL                  = EL(:);

%% Detect NaNs
indx                = find(isnan(EL));

%% Elevation 2 Phi
% Revert Elevations to standard elevation, i.e. reinstate front-back
% confusion!
selThetaRear        = EL>90 | EL<-90;
EL(selThetaRear)    = sign(EL(selThetaRear)).*(180-abs(EL(selThetaRear)));

%% Check for impossible double-polar coordinates
sel=abs(AZ)+abs(EL)>90;
if any(sel)
    MsgID = 'AudToolbox:azel2fart:NonExistentDoublePolarCoordinates';
    Msg = 'Non-existent Double Polar Coordinates! ABS(Azimuth)+ABS(Elevation) should be < 90';
    warning(MsgID,Msg);
end

%% Elevation 2 Phi
% Is front (mod 5 deg) or rear (mod 2.5 deg) speaker nearest
selrear         = mod(EL/2.5,2)>=0.5 & mod(EL/2.5,2)<=1.5;
% Convert to Speaker-number
Phi             = (EL+55+2.5.*selrear)./5+1+100.*selrear;
% Choose nearest speaker
Phi             = round(Phi);

%% Azimuth 2 Theta
% Theta               = round(real(asind(sind(AZ)./cosd(EL))));
Theta               = (real(asind(sind(AZ)./cosd(EL)))); % Removed "round" Jan 2009

%% Rear Position
% 
sel                 = selrear & ~selThetaRear;
Theta(sel)          = sign(Theta(sel)).*(180-abs(Theta(sel)));
sel                 = selrear & ~selThetaRear & AZ == 0;
Theta(sel)          = 180;

sel                 = ~selrear & selThetaRear;
Theta(sel)          = sign(Theta(sel)).*(180-abs(Theta(sel)));
sel                 = ~selrear & selThetaRear & AZ == 0;
Theta(sel)          = 180;

% Convert Theta of Back-Speakers
Theta(selrear) = -Theta(selrear);

%% Back to NaNs
Phi(indx)           = NaN;
Theta(indx)         = NaN;

%% Check for non-existent speakernumbers
realspkr        = [1:29 101:129 31 32];
sel             = ~ismember(Phi,realspkr);
if any(sel)
    msg = 'Invalid value for elevation:';
    msgid = 'FART:nonexistentElevation';
    warning(msgid,msg);
    disp([EL(sel) Phi(sel)])
end
Theta(sel)      = NaN;
Phi(sel)        = NaN;

%% Convert back to Matrix
Theta           = reshape(Theta,M,N);
Phi             = reshape(Phi,M,N);
