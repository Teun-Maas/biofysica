function Dat=loadtdt(tdtfile,nchan,nsample,ntrials)
% Load Data from Dat-files
%
% function [Dat]=loadDat(tdtfile,nchan,nsample,<ntrials>)
%
%   Reads Data from DatAFILE, based on nchan and nsample
%
%
%       tdtfile =   Datafile
%                       (e.g. XX0101<.Dat>)
%       nchan =     Number of channels
%                       (e.g. HEAD x, HEAD y, EYE x & EYE y, gives
%                       nchan = 4;)
%       nsample =   Number of samples
%                       nsample = trail duration (sec) * samplerate
%                       (e.g. trail duration is 2 sec. and the samplerate
%                       is 1000Hz, than nsample = 2000;)
%       ntrials =   Number of trials. If given loadDat checks if number
%                   of bytes corresponds to this number of trials.
%
%   See also: READCSV
%
% tomg Oct 2006
% Modified by: MarcW 2007

%% Initialization
nin             = 4;
ext             = '.tdt';
tdtfile         = fcheckext(tdtfile,ext);
tdtfile         = fcheckexist(tdtfile);
nbytesperfloat  = 4;

%% Check for number of trials
fid             = fopen(tdtfile,'r');
if ispc
    [Dat,n]     = fread(fid,[nchan*nsample,inf],'ushort');
else
    Dat         = fread(fid,[nchan*nsample,inf],'short','l');
end
ntrials = numel(Dat)/(nsample*nchan)
Dat             = reshape(Dat,nsample,nchan,ntrials);
Dat             = permute(Dat,[1 3 2]);
fclose(fid);

