function [B, kshift]=pa_spk_strfshiftmax(strf)
% [B,SHIFT] = PA_SPK_STRFSHIFTMAX(STRF)
%
% Shift the peak in the STRF to a central frequency.
%
%  Assumptions: 
% - STRF with 2.5 octave spectral width
% - 10 frequencies
%
% See also PA_SPK_RIPPLE2STRF, PA_SPK_PREDICT

% 2012 Marc van Wanrooij 
% e-mail: marcvanwanrooij@neural-code.com
% Original by: John van Opstal

r			= rms(strf,2);
[~,indx]	= maxpair(r);
nfreq		= size(strf,1);
kshift		= floor(nfreq/2)-indx;
B			= pa_spk_shiftmatc(strf,kshift);

function [mx,km] = maxpair(x)
n		= length(x);
y		= NaN(n-1,1);
for ii = 1:n-1,
    y(ii) = x(ii)+x(ii+1);
end
y(n)	= x(1)+x(n);
[mx,km] = max(y);

