function Signal = equalizer(x, filtname)
%  PA_EQUALIZER(X, FILTFILE)
%
%    Equalizing filtering of time sequence.
%
%    SIGNAL - Time sequence
%    FILTFILE - File name of the mat-file containing the filter (default:
%    visaton_equalizer.mat)
%
% See also FIR2, FILTFILT

% 2013 Marc van Wanrooij

%% Initialization
if nargin<2
filtname = which('hoop_visaton_equalizer.mat'); % general/average measured for visaton speakers in HOOP lab
end
load(filtname);
Signal      = filtfilt (b, 1, x);
