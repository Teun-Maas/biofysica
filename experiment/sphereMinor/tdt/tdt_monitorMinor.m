function handles = tdt_monitorMinor(handles)
% TDT_MONITOR
%
% Monitors TDT components
%
% This script can be called to create a figure window for TDT Active X Controls and to maintain it

set(handles.checkbox_RZ6connect,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTConnect));
set(handles.checkbox_RZ6load,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTLoad));
set(handles.checkbox_RZ6run,'value',bitget(handles.cfg.RZ6_1Status,handles.cfg.statTDTRun));

set(handles.text_RZ6,'String',['RZ6 #1: ' handles.cfg.sRZ6_1circuit '.rcx']);
