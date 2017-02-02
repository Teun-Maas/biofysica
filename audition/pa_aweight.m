function [W,A] = pa_aweight(F)
% [W,A] = PA_AWEIGHT(F)
%
% A-weighting function for sound levels at frequency F
%
%

% (c) 24 May 2011 Marc van Wanrooij

A		= getaweight(F);
A1000	= getaweight(1000);
W		= 20.*log10(A./A1000);

function A = getaweight(F)
%% Formula: http://en.wikipedia.org/wiki/DB(A)
A		= (12200^2 .* F.^4) ./ ((F.^2 +20.6^2) .*(F.^2 +12200^2).* (F.^2 +107.7^2).^0.5 .*(F.^2 +737.9^2).^0.5);
