function writeled(fid,LED,X,Y,Int,EventOn,Onset,EventOff,Offset)
% WRITELED(FID,LED,X,Y,INT,EVENTON,ONSET,EVENTOFF,OFFSET)
%
% Write a LED-stimulus line in an exp-file with file identifier FID.
%
% LED		- 'LED' for Hoop-LEDS or 'SKY' for Sky-LEDS
% X			- LED theta angle
% Y			- LED phi number (1-29 and 101-129)
% INT		- LED Intensity (0-255)
% EVENTON	- The Event that triggers the onset of the LED (0 - start of
% trial)
% ONSET		- The Time after the On Event (msec)
% EVENTOFF	- The Event that triggers the offset of the LED (0 - start of
% trial)
% OFFSET	- The Time after the Off Event (msec)
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

% (c) Marc van Wanrooij 2009
fprintf(fid,'%s\t%d\t%d\t \t%d\t%d\t%d\t%d\t%d\n',LED,X,Y,Int,EventOn,Onset,EventOff,Offset);
