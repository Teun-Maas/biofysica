function nirs = pa_nirs_downsample(nirs,varargin)
% NIRS = PA_NIRS_DOWNSAMPLE(NIRS,FDOWN);
%
% Downsample data in NIRS structure to FDOWN (default: 10 Hz).
%
% or
%
% NIRS = PA_NIRS_DOWNSAMPLE(O,FS,FDOWN);
%
% Downsample data in matrix O with sampling rate FS to FDOWN (default: 10 Hz).
%
%
%
% See also PA_NIRS_READ, PA_NIRS_PREPROCESS

% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

%% Initialization

if isstruct(nirs)
	O		= nirs.trial;
	M		= nirs.trial; % the data in a large matrix
	fsample = nirs.fsample;
	if nargin<2
		Fsdown = 10;
	elseif nargin==2
		Fsdown = varargin{1};
	end
else
	O = nirs;
	M = nirs;
	fsample = varargin{1};
	if nargin<3
		Fsdown = 10;
	elseif nargin==3
		Fsdown = varargin{2};
	end
end

%% downsample
[nchan,nsample] = size(M);

R				= fsample/Fsdown; % usually sample frequency = 250 Hz, so this resampling rate = 25
nrsample		= ceil(nsample/R);
Mdown			= NaN(nchan,nrsample);
for idxChan		= 1:nchan
	Mdown(idxChan,:) = decimate(M(idxChan,:),R);
% 		Mdown(idxChan,:) = downsample(M(idxChan,:),R);

end
M				= Mdown; % downsampled data matrix


%% Output
if isstruct(nirs)
nirs.trialdown		= M;
nirs.fsampledown	= Fsdown;
else
	nirs = M;
end