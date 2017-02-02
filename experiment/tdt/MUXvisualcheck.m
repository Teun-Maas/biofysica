function MUXvisualcheck
% MUXVISUALCHECK
%
% turn on all bits of thye PM2relay multiplexes one by one
%
% See also MUXSET

close all

%% TDT intialization
tdt_globals;		% TDT defaults, in cfg structure, includes circuit names for [RA16_1,RA16_2,RX6,RP2_1,RP2_2]
tdt_init;

%% Turn off everything
for muxIdx = 1:4
	MUX(RP2_1,muxIdx,0)
	MUX(RP2_2,muxIdx,0)
end

pause(1)

%% Disco
for muxIdx = 1:4
	for chnIdx = 1:16
		MUX(RP2_1,muxIdx,chnIdx);
		pause(.2)
	end
	MUX(RP2_1,muxIdx,0)	
end

for muxIdx = 1:4
	for chnIdx = 1:16
		MUX(RP2_2,muxIdx,chnIdx);
		pause(.2)
	end
	MUX(RP2_2,muxIdx,0)	
end

