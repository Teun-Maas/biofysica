function [W,A] = aweight(F)
% [W,A] = AWEIGHT(F)
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
% A		= (12200^2 .* F.^4) ./ ((F.^2 +20.6^2) .*(F.^2 +12200^2).* (F.^2 +107.7^2).^0.5 .*(F.^2 +737.9^2).^0.5);

%% From https://nl.mathworks.com/matlabcentral/fileexchange/46819-a-weighting-filter-with-matlab-implementation
c1 = 3.5041384e16;
c2 = 20.598997^2;
c3 = 107.65265^2;
c4 = 737.86223^2;
c5 = 12194.217^2;

f		= F.^2;
num		= c1*f.^4;
den		= ((c2+f).^2) .* (c3+f) .* (c4+f) .* ((c5+f).^2);
A		= num./den;
A		= A(:)';
