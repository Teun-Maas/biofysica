% signalIn = readWavFunc([name,] par)
%
% Read wav data from file. 
%
% INPUT:
%   name - name of wav-file (with path, generated with fullfile) [optional]
%   par  - parameter object / struct
%
% FIELDS FOR PAR:
%   sp_fs - desired sampling frequency (Hz)
%   wavFile - name of wav-file (used of readWavFunc is called with a single arguement)
%
% OUTPUT:
%   signalIn - samples of wav-file

% Change log:
%   Apr 2012, M. Milczynski - created
%   14 Jan 2013, PH - wav file can now be specified either by an input
%                   DataUnit (previous behavior), or by property 'wavFile' 
%                   (new behavior); might deprecate old behavior for
%                   release version
%   09 Apr 2013, PH - improved backwards compatibility with Matlab
function signalIn = readWavFunc(varargin)

if nargin == 1  % only 'par' specified
    par = varargin{1};
    name = par.wavFile;
elseif nargin == 2 % 'name' and 'par' specified
    name = varargin{1};
    par = varargin{2};
else
    error('Illegal number of arguments (required: 1 or 2).');
end

requiredFields = {'sp_fs'};
checkParamFields(par, requiredFields);
%name=strrep(name,'/','\');%%%%  change the backward slashes to forward slashes 
[signalIn, srcFs] = audioread(name);

signalIn = signalIn(:,1);

if srcFs ~= par.sp_fs
    signalIn = resample(signalIn, par.sp_fs, srcFs); 
end

