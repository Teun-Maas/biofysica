function Lp = pa_pressure2level(P)
% LP = PA_PRESSURE2LEVEL(P)
%
% Convert sound pressure P (Pa) to sound level LP (dB SPL).
%
%
% See also PA_OCT2BW, PA_LEVEl2PRESSURE

% (c) 24 May 2011 Marc van Wanrooij

Pref	= 20*10e-6; % 20 muPa RMS, usually considered human threshold of hearing at 1 kHz
Lp		= 20*log10(P./Pref);

function [W,A] = Aweight(F)
A = (12200^2 .* F.^4) ./ ((F.^2 +20.6^2) .*(F.^2 +12200^2).* (F.^2 +107.7^2).^0.5 .*(F.^2 +737.9^2).^0.5);

f = 1000;
A1000 = (12200^2 * f^4) / ((f^2 +20.6^2) *(f^2 +12200^2)* (f^2 +107.7^2)^0.5 *(f^2 +737.9^2)^0.5);

W = 20.*log10(A./A1000);