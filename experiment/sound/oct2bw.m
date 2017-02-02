
function F = oct2bw(F1,oct)
% F2 = OCT2BW(F1,OCT)
%
% Determine frequency F2 that lies OCT octaves from frequency F1
%

% (c) 2011-05-06 Marc van Wanrooij
F = F1 .* 2.^oct;
