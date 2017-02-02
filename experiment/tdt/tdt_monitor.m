% TDT_MONITOR
%
% Monitors TDT components
%
% This script can be called to create a figure window for TDT Active X Controls and to maintain it

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com

%% window
HF = findobj('Tag','ActXWin');
if isempty(HF)
	HF					= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
end
figure(HF);

% radio button handles & positions
pos			= [(0:2)'*50 ones(3,1) ones(3,1)*50 ones(3,1)*20];
HuiRA16_1	= nan(1,3);
HuiRA16_2	= nan(1,3);
HuiRX6		= nan(1,3);
HuiRP2_1	= nan(1,3);
HuiRP2_2	= nan(1,3);
for I_bit	= [cfg.statTDTConnect cfg.statTDTLoad cfg.statTDTRun]
	HuiRA16_1(I_bit) = uicontrol('Tag',['RA16_1Stat' num2str(I_bit)],'position',pos(I_bit,:)+[10 41 0 0],'style','radiobutton','enable','off');
	HuiRP2_1(I_bit)  = uicontrol('Tag',['RP2_1Stat' num2str(I_bit)],'position',pos(I_bit,:)+[10 21 0 0],'style','radiobutton','enable','off');
	HuiRP2_2(I_bit)  = uicontrol('Tag',['RP2_2Stat' num2str(I_bit)],'position',pos(I_bit,:)+[10 1 0 0],'style','radiobutton','enable','off');
end

% labels
str		= [{'Connect'};{'Load'};{'Run'}];
for I_bit	= [cfg.statTDTConnect cfg.statTDTLoad cfg.statTDTRun]
	uicontrol('position',pos(I_bit,:)+[0 71 0 0],'style','text','enable','on','string',str(I_bit));
end
uicontrol('position',[131 41 80 20],'style','text','enable','inactive','string','RA16 #1: ','HorizontalAlignment','right');
uicontrol('position',[131 21 80 20],'style','text','enable','inactive','string','RP2 #1: ','HorizontalAlignment','right');
uicontrol('position',[131 1 80 20],'style','text','enable','inactive','string','RP2 #2: ','HorizontalAlignment','right');
uicontrol('position',[211 41 200 20],'style','edit','enable','inactive','string',cfg.sRA16_1circuit);
uicontrol('position',[211 21 200 20],'style','edit','enable','inactive','string',cfg.sRP2_1circuit);
uicontrol('position',[211 1 200 20],'style','edit','enable','inactive','string',cfg.sRP2_2circuit);


%% set radio buttons according to status
for I_bit	= [cfg.statTDTConnect cfg.statTDTLoad cfg.statTDTRun]
	set(HuiRA16_1(I_bit),'value',bitget(RA16_1Status,I_bit),'enable','inactive');
	set(HuiRP2_1(I_bit),'value',bitget(RP2_1Status,I_bit),'enable','inactive');
	set(HuiRP2_2(I_bit),'value',bitget(RP2_2Status,I_bit),'enable','inactive');
end
drawnow