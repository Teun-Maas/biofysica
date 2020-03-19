function levelSPL = dBHL2SPL(frequency,levelHL)
% levelSPL = dBHL2SPL(frequency,levelHL) Convert a pure tone sound level
% from dB HL to dB SPL. Based on ANSI standard S3.6
% The inputs frequency and levelSPL must be either scalar or vectors of the
% same size. 
% This is the table used:
%     Hz		dB SPL 	dB HL
%     125 	45.0 	0
%     250 	27.0 	0
%     500 	13.5 	0
%     750	9.0 	0
%     1000	7.5		0
%     1500	7.5		0
%     2000	9.0		0
%     3000	11.5	0
%     4000	12.0	0
%     6000	16.0	0
%     8000	15.5	0
%	  See also: dBSPL2HL
% Written by: J.B. van der Heijdt
% Latest update: 2020-03-14

fTest = [125 250 500 750 1000 1500 2000 3000 4000 6000 8000];
if ~any(		fTest == frequency	)
	error('No conversion value defined for this value, please insert appropriate frequency');
end
for iFreq = 1:length(frequency) % in case of vector
	switch frequency(iFreq)
		case 125
			levelSPL(iFreq) = levelHL(iFreq) + 45.0;
		case 250
			levelSPL(iFreq) = levelHL(iFreq) + 27.0;
		case 500
			levelSPL(iFreq) = levelHL(iFreq) + 13.5;
		case 750 
			levelSPL(iFreq) = levelHL(iFreq) + 9.0;
		case 1000
			levelSPL(iFreq) = levelHL(iFreq) + 7.5;
		case 1500
			levelSPL(iFreq) = levelHL(iFreq) + 7.5;
		case 2000
			levelSPL(iFreq) = levelHL(iFreq) + 9.0;
		case 3000
			levelSPL(iFreq) = levelHL(iFreq) + 11.5;
		case 4000
			levelSPL(iFreq) = levelHL(iFreq) + 12.0;
		case 6000
			levelSPL(iFreq) = levelHL(iFreq) + 16.0;
		case 8000
			levelSPL(iFreq) = levelHL(iFreq) + 15.5;
	end
end

