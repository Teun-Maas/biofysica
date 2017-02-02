function [snd,Fs] = gennb(Dur,F0, Fs, varargin)
% Generate 1/3 octave-wide narrow-band 'Gaussian White Noise' Stimulus
%
% NB = GENNB (Dur,  F0, Fs)
%
% Generate narrow-band noise, with
% Dur       - duration of sound (s)
% F0        - central Frequency
% Fs        - sampling rate
%
% See also GENGWN
% http://nl.mathworks.com/help/audio/examples/octave-band-and-fractional-octave-band-filters.html

% 2016 Marc van Wanrooij

%% Initialization

if nargin<3
    Fs          = 48828.125; % TDT Nyquist sampling frequency (Hz)
end    
if nargin<2
    F0          = 2000; % Hz
end    

if nargin<1
	Dur = 0.15;
end
N    = round(Dur*Fs); % samples


%% Optional arguments
dspFlag       = keyval('display',varargin);
if isempty(dspFlag)
	dspFlag	= 1;
end
plee       = keyval('play',varargin);
if isempty(plee)
	plee	= 'n';
end

%% Create and Modify Signal
N				= round(N);
snd             = randn([N,1]);

%% Filter
BandsPerOctave	= 3;
N				= 8;   % Filter Order
f				= fdesign.octave(BandsPerOctave,'Class 1','N,F0',N,F0,Fs);
df				= design(f,'butter');
snd				= filter (df(1), snd);

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
if strcmpi(plee,'y');
	sndplay = envelope(snd',round(10*Fs/1000));
	p		= audioplayer(sndplay,Fs);
	playblocking(p);
end



