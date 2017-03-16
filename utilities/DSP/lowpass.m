function xfilt = lowpass(x,varargin)
%  XFILT = LOWPASS(X,NAME,VALUE)
%
%    Lowpass filtering of time sequence.
% 
%    X - Time sequence
%    'Fc'     - Cutoff frequency (default 3 kHz)
%    'Fs'     - sampling rate (default 48828 Hz)
%    'order' - Order of filter (default 100)
%
% See also FIR1,DESIGNFILTER, FILTFILT, HIGHPASS, GETPOWER

% (c) 2015 Marc van Wanrooij

%% Initialization
Fc = keyval('Fc',varargin);
if isempty(Fc)
	    Fc      = 500; % Cut-off frequency (Hz)
end
Fs = keyval('Fs',varargin);
if isempty(Fs)
	    Fs      = 48828; % Sample frequency
end
Fn = keyval('Fn',varargin);
if isempty(Fn)
	    Fn      = Fs/2;  % Nyquist frequency (Hz)
end
order = keyval('order',varargin);
if isempty(order)
    order  = 100;  % Filter order
end
if nargin<1
	x = randn(1000,1);
end
%% Filter
% f           = Fc/Fn;
% b           = fir1(order,f);
% xfilt      = filtfilt (b, 1, x);

lpFilt = designfilt('lowpassfir', 'FilterOrder', order, 'CutoffFrequency', ...
                    Fc, 'SampleRate', Fs, 'Window', 'hamming');
xfilt  = filtfilt (lpFilt, x);
% to compare filers: see fvtool
