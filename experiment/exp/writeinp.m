function writeinp(fid,chan)
% WRITEACQ(FID,EVENTON,ONSET)
%
% Write Acquisition line in an exp-file with file identifier FID.
%
% EVENTON	- The Event that triggers the acquisition (0 - start of
% trial)
% ONSET		- The Time after the On Event (msec)
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

% (c) Marc van Wanrooij 2009
if chan==1
fprintf(fid,'%s\t\t\t\t\t\t\t\t\n','INP1');
elseif chan==2
fprintf(fid,'%s\t\t\t\t\t\t\t\t\n','INP2');
end