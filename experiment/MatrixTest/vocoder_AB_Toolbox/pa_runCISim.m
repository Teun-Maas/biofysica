function Y = pa_runCISim(X,Fs,varargin)
% runCISim      Simplified interface to vocoderCISim.
%
% Usage:
%   OUT = runCISim(SIG,FS {,'paramName',paramVal...});
%
% Perform the following operations:
%   (1) resample to 17400
%   (2) call vocoderCISim
%   (3) resample back to FS
%   (4) normalize RMS to that of original signal
%   (5) ensure that size is the same.
%   (6) Add 3 dB of rove
%
% For additional processing options, see help on vocoderCISim.

nChans       = keyval('nChans',varargin); % Number of CI channels
nChans	     = 6;
%% Rove
dBRoveMax	= 3;
Fsnew		= 17400;
depth		= 2^15; % sound depth
X			= X/max(abs(X))*depth; % MW 'fix' fix
X			= X/depth;
XR			= resample(X,Fsnew,Fs); % resample the data

%% Add a little noise % MW: Why?
XR		= XR + ( (rand(size(XR))-0.5)*1e-25 );

%% Pre-emphasis % MW: What does this do?
A		= [1.0000   -1.5299    0.5453];
B		= [0.7688   -1.5376    0.7688];
XR		= filter(B,A, XR);

%% Vocde / Process
YR		= pa_vocoderCISim(XR,'nChans', nChans); 

%% De-emphasis
Binv	= A;
Ainv	= [0.7844   -1.5374    0.7533];
YR		= filter(Binv,Ainv, YR);

%% Resample to original
Y		= resample(YR,Fs,17400); 
if(length(Y) > length(X))
    Y = Y(1:length(X))';
else
    Y = [Y(:)' zeros(1,length(X)-length(Y))];
end

%% RMS adjustment
% to keep original sound level
rmsOrig		= rms(X);
rmsOut		= rms(Y);
Y			= Y*rmsOrig/rmsOut;

%% Roving
% This is just an experimental approach to induce uncertainty
% dbRoveLevel = rand(1,1)*dBRoveMax*2-dBRoveMax;
% Y			= Y*10^(dbRoveLevel/20) ;

%% Filter below 50 Hz to remove click
[b,a]	= butter(4,50/Fs*2,'high');
Y		= filter(b,a,Y);
%clf;plot(OUT)