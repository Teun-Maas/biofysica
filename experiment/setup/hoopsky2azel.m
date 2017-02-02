function [AZ,EL] = hoopsky2azel(X,Y)
% SKY2AZEL Transfrom SKY to double polar coordinates.
% [AZ,EL] = SKY2AZEL(R,Phi) transforms corresponding elements of data
%     stored in SKY coordinates R (Ring Number), Phi (Spoke Number)
%    to double polar coordinates (azimuth AZ and
%     elevation EL).  The arrays R and Phi must be the same size (or either
%     can be scalar). AZ and EL are returned in deg.
%
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
[M,N]       = size(X);
R           = X(:);
Phi         = Y(:);
sel = R==0;
if any(sel)
    MsgID   = 'AudToolbox:sky2azel:Ring0';
    Msg     = 'There is no Ring defined for center sky LED';
    warning(MsgID,Msg);
    R(sel)      = 1;
end
%% Ring
Dist        = [14 25.5 40 58.2 81.5 112 149.5];
R           = atand(Dist(R)/160)';

%% Center Fixation
R(Phi==0)   = 0;

%% Spoke
Phi         = -30*(Phi-1)+60;
sel         = Phi<0;
Phi(sel)    = Phi(sel)+360;

%% DoublePolar
[HV]        = pa_rphi2azel([R,Phi]);
AZ          = HV(:,1);
EL          = HV(:,2);

%% Convert back to MTX
AZ          = reshape(AZ,M,N);
EL          = reshape(EL,M,N);
