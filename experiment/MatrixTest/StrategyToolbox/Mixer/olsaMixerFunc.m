% [wavOut clip] = olsaMixerUnit(wav_1, wav_2, par)
%
% INPUT:
%   wav_i - vector containing wav data for input channel i
%   par   - paramteter struct / object
%
% FIELDS FOR PAR:
%   sp_fs    - input sampling rate
%   sensIn   - input sensitivity (dB SPL equivalent to sinusoid at digital FS)
%   lvlType  - 'rms' / 'peak' / 'rel'; 
%   lvlDb    - 2-element vector of levels  in dB
%   delays   - vector of onset delays for each input [s] 
%   wrap     - repeat shorter inputs to match duration of the longest
%              input? [1/0]
%   durFade  -  duration of cross-fade when wrapping signals [s]
%   clip     - read-only indicator: did clipping occur during the most recent
%              call of run()? [0/1]
%   nomSpeechRms - nominal speech input RMS on linear scale; compute actual RMS if empty
%   noiseSeek - initial seek position within noise [s];  -1 for random
%
% OUTPUT:
%    wavOut  - mixed signal (column vector)
%    clip    - clipping indicator [0/1]
%
% Change log:
%   04 Mar 2013, P.Hehrmann - created
%   09 Apr 2013, PH - improved backwards compatibility w. Matlab < 2011b
%   06 Aug 2013, PH - added: initial seek of noise; nominal speech RMS now optional
function [wavOut, clip] = olsaMixerFunc(varargin)

par = varargin{3};

wav = varargin([1,2]);

if isempty(wav{2}) 
    nWav = 1;
else
    nWav = 2; % speech and noise?
end

% make all wavs column vectors
for iWav = 1:nWav
    if size(wav{iWav},1) == 1 
        wav{iWav} = wav{iWav}(:);
    end
end

assert(all(cellfun(@(X__) isnumeric(X__) & isvector(X__), wav(1:nWav))), 'All wav inputs must be numerical vectors');
assert(length(par.lvlDb) == 2, 'Length of par.lvlDb must be 2.');
assert(isempty(par.delays) || (length(par.delays) == 2), 'Length of par.delays must be 0 or 2');
assert(ischar(par.lvlType), 'par.lvlType must be a string.');

% compute onset delay for each audio input
if isempty(par.delays)
    delays = zeros(1,nWav);
else
    assert(all(par.delays >= 0), 'Elements of par.delays must be non-negative.');
    delays = round(par.delays * par.sp_fs);
end

if par.noiseSeek >= 0
    nSeek = par.noiseSeek * par.sp_fs;
elseif par.noiseSeek == -1
    nSeek = round(rand*length(wav{2}));
end
delays(1) = delays(1)+nSeek;

% get level type for each input
lvlType = {'rms', par.lvlType};

% determine input signal lengths
lenWavIn = cellfun(@length, wav(1:nWav));
lenWavDelayed = lenWavIn + delays *2; % input length including delays && ADD HERE OFFSET DELAY!!!!!!!!

lenOut = lenWavDelayed(1);  % length of speech signal (wav_1)

% length of cross-fade in samples, and fade-in/out envelopes
lenFade = ceil(par.durFade * par.sp_fs);
envFadeOut = 0.5*cos(linspace(0,pi,lenFade))' + 0.5;
envFadeIn = 1-envFadeOut;

% determine input levels (prior to padding/wrapping)
lvlWav = NaN(1,2);

if isempty(par.nomSpeechRms)
    lvlWav(1) = 10*log10(mean(wav{1}.^2)) + par.sensIn;
else
    lvlWav(1) = 20*log10(par.nomSpeechRms) + par.sensIn; % <--
end

switch lower(lvlType{2})
    case 'rms'
        lvlWav(2) = 10*log10(mean(wav{2}.^2)) + par.sensIn;
    case 'peak'
        lvlWav(2) = 20*log10(max(abs(wav{2})))  + par.sensIn;
    case 'rel'
        lvlWav(2) = 0;
    otherwise
        error('Unknown level scaling type ''%s''', lvlType{2});
end

lvlWav(1) = 20*log10(rms(wav{1})/par.nomSpeechRms) + par.sensIn;
lvlWav(2) = 20*log10(rms(wav{2})/par.nomSpeechRms) + par.sensIn;
%[par.sensIn par.nomSpeechRms rms(wav{1}) rms(wav{2})]

% find wavs that need to be wrapped / padded
needsLengthening = (lenWavDelayed < lenOut);

for iWav = 1:nWav
    if needsLengthening(iWav)
        if par.wrap % wrap signal
            nRep = ceil( (lenOut-delays(iWav))/(lenWavIn(iWav)-lenFade) - 1 );
            wavCross = envFadeOut .* wav{iWav}(end-lenFade+1:end) + envFadeIn .* wav{iWav}(1:lenFade);
            wav{iWav} = [zeros(delays(iWav),1); wav{iWav}(1:end-lenFade); ...
                repmat([wavCross; wav{iWav}(lenFade+1:end-lenFade)], nRep, 1 )];
            wav{iWav}(lenOut+1:end) = [];
        else % zero-pad signal
            wav{iWav} = [zeros(delays(iWav),1); wav{iWav}; zeros(lenOut-lenWavDelayed(iWav),1)];
        end
    else % truncate signal
        wav{iWav} = [zeros(delays(iWav),1); wav{iWav}(1:lenOut-delays(iWav)*2); zeros(delays(iWav),1)];
        wav{iWav}(end+1:lenOut) = [];
    end    
end

assert(all(cellfun(@(W__) length(W__) == lenOut, wav(1:nWav))), 'All wavs must have length lenMax by now.')

% max presentable sound = 81;
% add scaled inputs
wavOut = zeros(lenOut, 1);
for iWav = 1:nWav % 2x: speech and noise
    wavOut = wavOut + wav{iWav} * 10^((par.lvlDb(iWav)-lvlWav(iWav)) / 20);
%    [iWav par.lvlDb(iWav) lvlWav(iWav) par.lvlDb(iWav)-lvlWav(iWav) lvlWav(iWav)+par.lvlDb(iWav)-lvlWav(iWav) max(abs(wavOut(:)))]
end

% remove initial samples up to seek position
wavOut(1:nSeek) = [];

maxAbsOut = max(abs(wavOut));
clip = (maxAbsOut > 1);
if clip
    warning('Clipping occured. Maximum output amplitude %.2f (%.2fdB FS, equiv. %.2fdB SPL)', maxAbsOut, 20*log10(maxAbsOut), 20*log10(maxAbsOut)+par.sensIn);
end

end