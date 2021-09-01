% [wavOut clip] = audioMixerUnit(wav_1, ..., wav_n, par)
% Mix arbitrary number of audio inputs signals wav_i. For every signal, the
% target level as well as the target level type have to be specified 
% (abs. rms/abs. peak/rel. to input). par.sensIn defines the peak SPL 
% equivalent to 0dB FS. par.wrap controls the behaviour of input signals
% are of unequal length (warp around or zero-pad to match duration of the 
% primary signal). For par.wrap = 1, par.durFade controls the duration of 
% a cosine-shaped cross-fade between the end and beginning.
%
% INPUT:
%   wav_i - vector containing wav data for input channel i
%   par   - paramteter struct / object
%
% FIELDS FOR PAR:
%   sp_fs    - input sampling rate
%   sensIn   - input sensitivity: dB SPL corresponding to digital RMS amp. (0dB re.) 1
%                  [equivalently: dB peak corresponding to digital peak amp. 1]  
%   lvlType  - string or cell array of strings: 'rms' / 'peak' / 'rel'; 
%              if string, same type is used for all channels;
%              if cell aray, the size must match n
%   lvlDb    - n-vector of levels per channel in dB; for types 'rms' and
%              'peak', lvlDb(i) is in dB SPL; for type 'rel', lvlDb(i) is
%              in dB relative to the input level.
%   delays   - vector of onset delays for each input [s] 
%   primaryIn - which input determines the length of the mixed output
%               (1..nInputs, or []); if [], the longtest input (including 
%               delay) is chosen as primary.
%   wrap     - repeat shorter inputs to match duration of the primary
%              input? [1/0]
%   durFade  -  duration of cross-fade when wrapping signals [s]
%
% OUTPUT:
%    wavOut  - mixed signal (column vector)
%    clip    - clipping indicator [0/1]
%
% Change log:
%   29/08/2012, P.Hehrmann - created
%   31/08/2012, P.Hehrmann - bug fix: determine input levels before
%                            wrapping/padding
%   12/09/2012, P.Hehrmann - bug fix (unwanted error occurred in case of clipping)
%   09/12/2012, PH - added 'primaryIn' and 'delays' functionality
%   11/12/2012, PH - convenience fix: make all inputs column vectors
function [wavOut clip] = audioMixerFunc(varargin)

par = varargin{end};

nWav = nargin-1;
wav = varargin(1:nWav);
% make all wavs column vectors
for iWav = 1:nWav
    if isrow(wav{iWav})
        wav{iWav} = wav{iWav}(:);
    end
end

assert(all(cellfun(@(X__) isnumeric(X__) & isvector(X__), wav(1:nWav))), 'wav_1..wav_n must be numerical vectors');
assert(length(par.lvlDb) == nWav, 'Length of par.lvlDb must equal the number of audio inputs.');
assert(isempty(par.delays) || (length(par.delays) == nWav), 'Length of par.delays must 0 or equal the the number of audio inputs. ' )
assert(ischar(par.lvlType) || iscellstr(par.lvlType), 'par.lvlType must be a string of cell array of strings.')
assert(isempty(par.primaryIn) || (~mod(par.primaryIn,1) && (par.primaryIn <= nWav) && par.primaryIn > 0),...
       'par.primaryIn must be empty or an integer less or equal to the number of audio inputs.');

% compute onset delay for each audio input
if isempty(par.delays)
    delays = zeros(1,nWav);
else
    assert(all(par.delays >= 0), 'Elements of par.delays must be non-negative.');
    delays = round(par.delays * par.sp_fs);
end

% get level type for each input
if ischar(par.lvlType)
    lvlType = repmat({par.lvlType},1,nWav);
else
    lvlType = par.lvlType;
end

% determine input signal lengths
lenWavIn = cellfun(@length, wav(1:nWav));
lenWavDelayed = lenWavIn + delays; % input length including delays

% determine output length
if isempty(par.primaryIn)
    lenOut = max(lenWavDelayed);
else
    lenOut = lenWavDelayed(par.primaryIn);
end

% length of cross-fade in samples, and fade-in/out envelopes
lenFade = ceil(par.durFade * par.sp_fs);
envFadeOut = 0.5*cos(linspace(0,pi,lenFade))' + 0.5;
envFadeIn = 1-envFadeOut;

% determine input levels (prior to padding/wrapping)
lvlWav = NaN(1,nWav);
for iWav = 1:nWav
    switch lower(lvlType{iWav})
        case 'rms'
            lvlWav(iWav) = 10*log10(mean(wav{iWav}.^2)) + par.sensIn;
        case 'peak'
            lvlWav(iWav) = 20*log10(max(abs(wav{iWav})))  + par.sensIn;
        case 'rel'
            lvlWav(iWav) = 0;
        otherwise 
            error('Unknown level scaling type ''%s'' at index %d', lvlType{iWav}, iWav);
    end
end

% find wavs that need to be wrapped / padded
needsLengthening = (lenWavDelayed < lenOut);

for iWav = 1:nWav
    % RETHINK nRep with delays!
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
        wav{iWav} = [zeros(delays(iWav),1); wav{iWav}(1:lenOut-delays(iWav))];
        wav{iWav}(end+1:lenOut) = [];
    end
    
end

assert(all(cellfun(@(W__) length(W__) == lenOut, wav)), 'All wavs must have length lenMax by now.')

% add scaled inputs
wavOut = zeros(lenOut, 1);
for iWav = 1:nWav
    wavOut = wavOut + wav{iWav} * 10^((par.lvlDb(iWav)-lvlWav(iWav)) / 20);
end

maxAbsOut = max(abs(wavOut));
clip = (maxAbsOut > 1);
if clip
    warning('Clipping occured. Maximum output amplitude %.2f (%.2fdB FS, equiv. %.2fdB SPL)', maxAbsOut, 20*log10(maxAbsOut), 20*log10(maxAbsOut)+par.sensIn);
end

end