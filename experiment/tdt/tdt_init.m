% TDTINIT
%
% Initializes TDT system, checks and opens graphical monitor
%
% See also TDT_GLOBALS, TDT_MONITOR


%% Initialization
fprintf('-----------------------------------\n');
fprintf('Initializing TDT\n');
t = tic;

%% Active X Control/Objects
HF					= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
[zBus, err(1)]		= ZBUS(1); % zBus, number of racks
[RP2_1, err(2)]		= RP2(1,cfg.RP2_1circuit); % Real-time processor 1
[RP2_2, err(3)]		= RP2(2,cfg.RP2_2circuit); % Real-time processor 2
[RA16_1, err(4)]	= RA16(1,cfg.RA16_1circuit); % Real-time acquisition
[PA5_1, err(5)]		= PA5(1); % Programmable attenuator 1
[PA5_2, err(6)]		= PA5(2); %  Programmable attenuator 2
[PA5_3, err(7)]		= PA5(3); %  Programmable attenuator 3
[PA5_4, err(8)]		= PA5(4); %  Programmable attenuator 4
for muxIdx = 1:4
	MUX(RP2_1,muxIdx);
	MUX(RP2_2,muxIdx);
end

%% TDT status
RA16_1Status	= RA16_1.GetStatus;
RP2_1Status		= RP2_1.GetStatus;
RP2_2Status		= RP2_2.GetStatus;
tdt_monitor;
toc(t)
fprintf('-----------------------------------\n');

%% Configuration
cfg.RP2_1	= RP2_1;
cfg.RP2_2	= RP2_2;
cfg.RA16_1	= RA16_1;
cfg.zBus	= zBus;
cfg.PA5_1	= PA5_1;
cfg.PA5_2	= PA5_2;
cfg.PA5_3	= PA5_3;
cfg.PA5_4	= PA5_4;

