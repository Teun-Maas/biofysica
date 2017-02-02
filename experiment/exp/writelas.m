function writelas(fid,LAS,ID,EventOn,Onset,EventOff,Offset,Event)
% WRITELAS(FID,LAS,ID,EVENTON,ONSET,EVENTOFF,OFFSET)
%
% Write a LAS-stimulus line in an exp-file with file identifier FID.
%
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

% (c) Marc van Wanrooij 2009
fprintf(fid,'%s\t\t\t%d\t\t%d\t%d\t%d\t%d\t%d\n',LAS,ID,EventOn,Onset,EventOff,Offset,Event);
