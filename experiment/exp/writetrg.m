function writetrg(fid,Edge,Bit,EventOn,Onset,ID)
% WRITELED(FID,EDGE,BIT,EVENTON,ONSET,ID)
%
% Write a trigger-event line in an exp-file with file identifier FID.
%
% EDGE		- EdgeFlag: 1 - Button Press, 2 - Button Release
% BIT		- MicroController Bit
% EVENTON	- The Event that triggers the onset of the TRG (0 - start of
% trial)
% ONSET		- The Time after the On Event (msec)
% ID		- The ID of the TRG-event
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

% (c) Marc van Wanrooij 2009
fprintf(fid,'%s\t\t\t%d\t%d\t%d\t%d\t\t\t%d\n','TRG0',Edge,Bit,EventOn,Onset,ID);