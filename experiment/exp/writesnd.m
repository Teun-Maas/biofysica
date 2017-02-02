function writesnd(fid,SND,X,Y,ID,Int,EventOn,Onset,Extra)
% WRITESND(FID,SND,X,Y,ID,INT,EVENTON,ONSET)
%
% Write a SND-stimulus line in an exp-file with file identifier FID.
%
% SND		- 'SND1' or 'SND2'
% X			- SND theta angle
% Y			- SND phi number (1-29 and 101-129)
% ID		- Sound ID/attribute (000-999)
% INT		- SND Intensity (0-100)
% EVENTON	- The Event that triggers the onset of the SND (0 - start of
% trial)
% ONSET		- The Time after the On Event (msec)
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

if nargin<9
fprintf(fid,'%s\t%d\t%d\t%d\t%.1f\t%d\t%d\n',SND,X,Y,ID,Int,EventOn,Onset);
else
fprintf(fid,'%s\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\n',SND,X,Y,ID,Int,EventOn,Onset,Extra);
end	