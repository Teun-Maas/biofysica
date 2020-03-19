function levelHL = dBSPL2HL(frequency,levelSPL)
% levelHL = dBSPL2HL(frequency,levelSPL) Convert a pure tone sound level
% from dB SPL to dB HL. Based on ANSI standard S3.6.
% The inputs frequency and levelSPL must be either scalar or vectors of the
% same size. 
% This is the table used:
%     Hz	dB SPL 	dB HL
%     125 	45.0 	0
%     250 	27.0 	0
%     500 	13.5 	0
%     750	9.0 	0
%     1000	7.5	0
%     1500	7.5	0
%     2000	9.0	0
%     3000	11.5	0
%     4000	12.0	0
%     6000	16.0	0
%     8000	15.5	0
%	 See also: dBHL2SPL
% Written by: J.B. van der Heijdt
% Latest update: 2019-12-17


fTest = [125 250 500 750 1000 1500 2000 3000 4000 6000 8000];
levelHL = nan; % preallocate

for iFreq = 1:length(frequency) % in case of vector
	% check inputs
	if ~any(		fTest == frequency(iFreq)	)
		warning('No conversion value defined for this value, please insert appropriate frequency');
		levelHL(iFreq) = nan;
	else
		switch frequency(iFreq)
			case 125
				levelHL(iFreq) = levelSPL(iFreq) - 45.0;
			case 250
				levelHL(iFreq) = levelSPL(iFreq) - 27.0;
			case 500
				levelHL(iFreq) = levelSPL(iFreq) - 13.5;
			case 750 
				levelHL(iFreq) = levelSPL(iFreq) - 9.0;
			case 1000
				levelHL(iFreq) = levelSPL(iFreq) - 7.5;
			case 1500
				levelHL(iFreq) = levelSPL(iFreq) - 7.5;
			case 2000
				levelHL(iFreq) = levelSPL(iFreq) - 9.0;
			case 3000
				levelHL(iFreq) = levelSPL(iFreq) - 11.5;
			case 4000
				levelHL(iFreq) = levelSPL(iFreq) - 12.0;
			case 6000
				levelHL(iFreq) = levelSPL(iFreq) - 16.0;
			case 8000
				levelHL(iFreq) = levelSPL(iFreq) - 15.5;
		end
	end
	
end

