function [snd,Fs] = gentone(Freq, dur, varargin)
% Generate TONE Stimulus
%
% STM = GENTONE (FREQ, DUR, NENV)
%
% Generate a sine-shaped tone, with
% dur       - duration                 [150]      ms
% Freq      - Frequency of tone                 [4000]      Hz
%
% PA_GENTONE(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
% 'nenvelope' - number of samples in envelope     [250]       samples
%             (head and tail are both 'NEnvelope')
% 'Fs'			- Sample frequency                  [48828.125] Hz
% 'display'		- set to 1 to display stuff         [0]
% 'phase'		- set phase							[0]
%
%
% For example:
%   Sine440 = gentone(440,150);
%   wavplay(Sine440,48828.125)
%   % Standard 150ms 440Hz tone
%
%   Sine3000 = gentone(3000,150,'nenvelope',7200/2,'Fs',50000,'display',1);   
%   wavplay(Sine3000,50000)
%   % 3kHz tone with onset and offset ramp with a sample rate of 50kHz. It 
%   % also produces an output screen of the generated tone in time and
%   % frequency domain.
%
% See also WRITEWAV, LOWPASSNOISE, HIGHPASSNOISE, GENGWN
%

%% Initialization
if nargin<1
    Freq = 4000; % Hz
end    
if nargin<2
    dur = 150; % samples
end
dspFlag		= keyval('display',varargin,0);
plee		= keyval('play',varargin,false);
Fs			= keyval('Fs',varargin,48828.125);
phi			= keyval('phase',varargin,0);
N			= round(dur/1000*Fs);
Fn			= Fs/2;

%% Create and Modify Signal
sig				= 0:N-1;
sig             = sig/Fs;
snd             = sin(2*pi*Freq*sig+phi);

%% Optional Graphics
if dspFlag
    figure;
    disp(['>> ' upper(mfilename) ' <<']);
    subplot(211)
    plot(snd)
    xlabel('Sample number')
    ylabel('Amplitude (a.u.)');
    
    Nfft = 2^10;
    Nnyq = Nfft/2;
    s = fft(snd,Nfft);
    s = abs(s);
    s = s(1:Nnyq);
    f = (0:(Nnyq-1))/(Nnyq)*Fn;
    subplot(212)
    loglog(f,s,'ko-','MarkerFaceColor','w');
    set(gca,'Xtick',[1 2 4 6 12 24]*1000,'XtickLabel',[1 2 4 6 12 24]);
    xlabel('Frequency (Hz');
    ylabel('Amplitude (au)');
end

if plee
	sndplay = envelope(snd',round(10*Fs/1000));
	p		= audioplayer(sndplay,Fs);
	playblocking(p);
end