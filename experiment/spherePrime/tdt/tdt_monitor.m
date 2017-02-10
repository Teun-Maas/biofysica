function handles = tdt_monitor(handles)
% TDT_MONITOR
%
% Monitors TDT components
%
% This script can be called to create a figure window for TDT Active X Controls and to maintain it

set(handles.checkbox_RA16connect,'value',bitget(handles.cfg.RA16_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RA16load,'value',bitget(handles.cfg.RA16_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RA16run,'value',bitget(handles.cfg.RA16_1Status,handles.cfg.statTDTRun));

set(handles.checkbox_RP21connect,'value',bitget(handles.cfg.RP2_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RP21load,'value',bitget(handles.cfg.RP2_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RP21run,'value',bitget(handles.cfg.RP2_1Status,handles.cfg.statTDTRun));

set(handles.checkbox_RP22connect,'value',bitget(handles.cfg.RP2_2Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RP22load,'value',bitget(handles.cfg.RP2_2Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RP22run,'value',bitget(handles.cfg.RP2_2Status,handles.cfg.statTDTRun));
set(handles.text_RA16,'String',['RA16 #1: ' handles.cfg.sRA16_1circuit '.rcx']);
set(handles.text_RP21,'String',['RP2 #1: ' handles.cfg.sRP2_1circuit '.rcx']);
set(handles.text_RP22,'String',['RP2 #2: ' handles.cfg.sRP2_2circuit '.rcx']);
