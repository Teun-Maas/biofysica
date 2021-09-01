function OUT = runCISim(SIG,FS,varargin)
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

% Rove
dBRoveMax = 3;

SIGR=resample(SIG/2^15,17400,FS);

%==============================================
% Add a little noise
SIGR = SIGR + ( (rand(size(SIGR))-0.5)*1e-25 );

%==============================================
% Pre-emphasis
A = [1.0000   -1.5299    0.5453];
B = [0.7688   -1.5376    0.7688];

SIGR = filter(B,A, SIGR);


%==============================================
% Process
OUTR=vocoderCISim(SIGR,varargin{:}); 

%==============================================
% De-emphasis
Binv = A;
Ainv = [0.7844   -1.5374    0.7533];

OUTR = filter(Binv,Ainv, OUTR);

% TODO: fix if length of resampled is less than original
OUT=resample(OUTR,FS,17400); 

if(length(OUT) > length(SIG))
    OUT=OUT(1:length(SIG))';
else
    OUT = [OUT(:)' zeros(1,length(SIG)-length(OUT))];
end

% RMS adjustment
rmsOrig = sqrt(sum(SIG.*SIG));
rmsOut = sqrt(sum(OUT.*OUT));

OUT = OUT * (rmsOrig/rmsOut);

% Roving
dbRoveLevel = rand(1,1)*dBRoveMax*2-dBRoveMax;
OUT = OUT * 10^(dbRoveLevel/20) ;

% Filter below 50 Hz to remove click
[b,a] = butter(4,50/FS*2,'high');
OUT = filter(b,a,OUT);
%clf;plot(OUT)