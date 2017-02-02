function xfilt = highpass(x,varargin)
%  XFILT = HIGHPASS(X,NAME,VALUE)
%
%    Highpass filtering of time sequence.
% 
%    X - Time sequence
%    'Fc'     - Cutoff frequency (default 3 kHz)
%    'Fs'     - sampling rate (default 48828 Hz)
%    'order' - Order of filter (default 100)
%
% See also FIR1,DESIGNFILTER, FILTFILT, LOWPASS, GETPOWER

% (c) 2015 Marc van Wanrooij

%% Initialization
Fc = keyval('Fc',varargin);
if isempty(Fc)
	    Fc      = 3000; % Cut-off frequency (Hz)
end
Fs = keyval('Fs',varargin);
if isempty(Fs)
	    Fs      = 48828; % Sample frequency
end
order = keyval('order',varargin);
if isempty(order)
    order  = 50;  % Filter order
end
if nargin<1
	x = randn(1000,1);
end
%% Filter
% Old
% f           = Fc/(Fs/2);
% b           = fir1(order,f,'high');
% xfilt      = filtfilt (b, 1, x);

% New
hpFilt = designfilt('highpassfir', 'FilterOrder', order, 'CutoffFrequency', ...
                    Fc, 'SampleRate', Fs, 'Window', 'hamming');
xfilt  = filtfilt (hpFilt, x);
% to compare filers: see fvtool
