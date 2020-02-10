function [snd,Fs] = gengwn(Dur, order, Fc, Fs, Fh, varargin)
% Generate Gaussian White Noise Stimulus
%
% GWN = GENGWN (Dur, order, Fc, Fs, Fh)
%
% Generate Gaussian White Noise Stimulus, with
% Dur       - duration of sound (s)
% Order     - order of filter
% Fc        - low-pass cut-off Frequency
% Fs        - sampling rate
% Fh		- high-pass cutoff
%
% For example:
%   Dur		= 0.15;
%   stm		= gengwn(Dur,100,20000,25000)
%   fname	= 'BB.wav';
%   writewav(stm,fname);
%
% will generate a broad-band noise between 500 (default) and 20 kHz with
% duration 150 msec (7500 samples / (25000 samples/sec*2) *1000 m). The
% on- and offset ramp each contain 250 samples = 250/50000 sec = 5 msec.
% This noise will be stored in the WAV-file 'BB.wav'.
%
% See also LOWPASS, HIGHPASS

% 2007 Marc van Wanrooij

%% Initialization
if nargin<5
    Fh          = 500;
end
if nargin<4
    Fs          = 48828.125; % TDT Nyquist sampling frequency (Hz)
end    
if nargin<3
    Fc          = 20000; % Hz
end    
if nargin<2
    order       = 100; % samples
end
if nargin<1
	Dur = 0.15;
end
N    = round(Dur*Fs); % samples
if order>N/6
	order = round(N/3)-2;
end

%% Optional arguments
dspFlag       = keyval('display',varargin);
if isempty(dspFlag)
	dspFlag	= 0;
end
plee       = keyval('play',varargin);
if isempty(plee)
	plee	= 'n';
end

%% Create and Modify Signal
N				= round(N);
snd             = randn([N,1]);
% Low-pass filter
snd             = lowpass(snd,'Fc',Fc,'Fs',Fs,'order', order);
% High-pass filter
snd             = highpass(snd,'Fc',Fh,'Fs',Fs,'order', order);

%% Optional Graphics
if dspFlag
    figure;
    disp('>> GENGWN <<');
    subplot(211)
    plot(snd)
    xlabel('Sample number')
    ylabel('Amplitude (a.u.)');
    
    subplot(212);
    getpower(snd,Fs,'display',1);
	xlim([50 22000])
end

%% Play
if strcmpi(plee,'y')
	sndplay = envelope(snd',round(10*Fs/1000));
	p		= audioplayer(sndplay,Fs);
	playblocking(p);
end