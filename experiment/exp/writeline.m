function writeline(fid,Type,X,Y,ID,Int,EventOn,On,EventOff,Off,Event)
% WRITETRL(FID,TRL)
%
% Write line in an exp-file with file identifier FID.
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\n',Type,X,Y,ID,Int,EventOn,On,EventOff,Off,Event);