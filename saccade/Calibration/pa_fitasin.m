function Par = pa_fitasin(X,Y,Inipar)
% PA_FITSIN   Fit an arcsine through data.
%
%   PAR = FITASIN(X,Y) returns parameters of the arcine a*asin(b*X+c)+d,
%   in the following order: [a b c d].
%
%   See also FMINSEARCH, NORM

% 2011, Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com


X           = X(:);
Y           = Y(:);
if nargin<3
  Inipar(1) = 15000;
  Inipar(2) = pi/180;
  Inipar(3) = 0;
  Inipar(4) = 0;
end
Par         = fminsearch(@asinerr,Inipar,[],X,Y);

function err =  asinerr(Par,X,Y)
%ASINERR   Determines error between experimental data and calculated arcsine.
%   ASINERR is used by the function FITASIN, and cannot be used by itself.
%
%   [ERR]=ASINERR(PAR,X,Y) returns the error between the calculated parameters
%   PAR, given by FITASIN and the parameters given by experimental data X and Y.
err = norm(Y-pa_asin(X,Par));


