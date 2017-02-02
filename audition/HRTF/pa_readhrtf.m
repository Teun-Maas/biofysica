function hrtf = pa_readhrtf(hrtffile)
% HRTF = PA_READHRTF(HRTFFILE)
%
% Read hrtf-file and extract Time Series data in HRTF
% HRTF is a (k x m x n)-matrix, where
%   k = number of samples
%   m = number of trials/locations
%   n = number of microphones/channels
%
% See also READCSV, READDAT

% Copyright 2008
% Marc van Wanrooij

%% Initialization
% Check Inputs
if nargin <1
    hrtffile        = '';
end
% Check Files
hrtffile            = pa_fcheckext(hrtffile,'.hrtf');
hrtffile            = pa_fcheckexist(hrtffile,'.hrtf');
csvfile             = pa_fcheckext(hrtffile,'.csv');

% Set Parameters/Constants
nbytesperfloat      = 4;

%% Read CSV-file
[~,~,mlog]     = pa_readcsv(csvfile);
% to determine number of trials

% nrepeats			= exp(3);
% ntrials             = exp(4);
% and number of channels
sel                 = ismember(mlog(:,5),[6 7]); % 6 and 7 are Inp1 and Inp2
nchan               = length(unique(mlog(sel,5)));

sel                 = ismember(mlog(:,5),6); % 6 and 7 are Inp1 and Inp2
ntrials				= sum(sel);
% mlog
%% Read HRTF-file
fid             = fopen(hrtffile,'r');
% First determine number of samples
fseek(fid,0,'eof');
nbytes          = ftell(fid);
nsample         = nbytes./(nchan*ntrials*nbytesperfloat);
% Get HRTF data
frewind(fid)
if ispc
    hrtf         = fread(fid,inf,'float');
else %MAC
    hrtf         = fread(fid,inf,'float','l');
end
% Some Reshaping

hrtf             = reshape(hrtf,nsample,nchan,ntrials);
% hrtf             = permute(hrtf,[1 3 2]); % Costs memory
% And close
fclose(fid);
