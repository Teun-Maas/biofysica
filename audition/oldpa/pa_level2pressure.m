function P = pa_level2pressure(Lp,BW)
% P = PA_LEVEL2PRESSURE(LP)
%
% Convert sound level LP (dB SPL) to sound pressure P (Pa).
%
% P = PA_LEVEL2PRESSURE(LP,BW)
%
% Convert sound level (dB SPL) to sound pressure per frequency (Pa/Hz) for
% a sound of bandwidth BW (default: sound pressure of entire sound).
%
% See also PA_OCT2BW, PA_PRESSURE2LEVEl

% (c) 24 May 2011 Marc van Wanrooij
if nargin<1
	Lp		= 60;
end
if nargin<2
% BW		= 20000;
	BW = 1;
end
Pref	= 20*10e-6; % 20 muPa RMS, usually considered human threshold of hearing at 1 kHz
P		= 10^(Lp/20)*Pref/BW;