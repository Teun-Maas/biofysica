function handles = tdt_initMinor(handles)
% TDTINIT
%
% Initializes TDT system, checks and opens graphical monitor
%
% See also TDT_GLOBALS, TDT_MONITOR



%% Active X Control/Objects
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
zBus		= ZBUS(1); % zBus, number of racks
RZ6_1		= RZ6(1,handles.cfg.RZ6_1circuit); % Real-time acquisition

Fs					= RZ6_1.GetSFreq;
handles.cfg.RZ6Fs	= Fs;
for muxIdx = 1:2
	MUX(RZ6_1,muxIdx);
end

%% TDT status
handles.cfg.RZ6_1Status	= RZ6_1.GetStatus;
handles					= tdt_monitorMinor(handles);

%% Configuration
handles.cfg.RZ6_1	= RZ6_1;
handles.cfg.zBus	= zBus;
