function handles = tdt_init(handles)
% TDTINIT
%
% Initializes TDT system, checks and opens graphical monitor
%
% See also TDT_GLOBALS, TDT_MONITOR


t = tic;

%% Active X Control/Objects
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
[zBus, err(1)]		= ZBUS(1); % zBus, number of racks
[RP2_1, err(2)]		= RP2(1,handles.cfg.RP2_1circuit); % Real-time processor 1
[RP2_2, err(3)]		= RP2(2,handles.cfg.RP2_2circuit); % Real-time processor 2
[RA16_1, err(4)]	= RA16(1,handles.cfg.RA16_1circuit); % Real-time acquisition
[PA5_1, err(5)]		= PA5(1); % Programmable attenuator 1
[PA5_2, err(6)]		= PA5(2); %  Programmable attenuator 2
[PA5_3, err(7)]		= PA5(3); %  Programmable attenuator 3
[PA5_4, err(8)]		= PA5(4); %  Programmable attenuator 4
for muxIdx = 1:4
	MUX(RP2_1,muxIdx);
	MUX(RP2_2,muxIdx);
end

%% TDT status
handles.cfg.RA16_1Status	= RA16_1.GetStatus;
handles.cfg.RP2_1Status		= RP2_1.GetStatus;
handles.cfg.RP2_2Status		= RP2_2.GetStatus;
handles					= tdt_monitor(handles);
toc(t)

%% Configuration
handles.cfg.RP2_1	= RP2_1;
handles.cfg.RP2_2	= RP2_2;
handles.cfg.RA16_1	= RA16_1;
handles.cfg.zBus	= zBus;
handles.cfg.PA5_1	= PA5_1;
handles.cfg.PA5_2	= PA5_2;
handles.cfg.PA5_3	= PA5_3;
handles.cfg.PA5_4	= PA5_4;
