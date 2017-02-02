function writeacq(fid,EventOn,Onset,Dur)
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
if nargin<4
fprintf(fid,'%s\t\t\t\t\t%d\t%d\n','ACQ',EventOn,Onset);
else
fprintf(fid,'%s\t\t\t\t\t%d\t%d\t%d\n','ACQ',EventOn,Onset,Dur);
end