function handles = datShow(handles)


%% load data
Hchan		= 3; Vchan		= 1; Fchan		= 2; % set channels in correct order
S			= load(handles.cfg.calfile,'-mat');
dat			=  handles.data(1).raw;
H			= dat(:,Hchan); H = H(:);
V			= dat(:,Vchan); V = V(:);
F			= dat(:,Fchan); F = F(:);
DAT			= [H V F]';

%% calibrate
[az,el]		= pa_calib(DAT,S); % calibrate: from Volts 2 Deg
az          = filtfilt(handles.cfg.lpFilt,az); % filter
el          = filtfilt(handles.cfg.lpFilt,el); % filter

t			= (1:numel(az))/handles.cfg.medusaFs; % time (s)
[vel,smv]	= getvel(az',el',handles.cfg.medusaFs,0.01); % get velocities

%% plot
if isfield(handles.cfg,'hcurdat')
	n = numel(handles.cfg.hcurdat);
	for figIdx = 1:n
		if ismember(figIdx,[1:3 5])
			set(handles.cfg.hcurdat(figIdx),'LineWidth',1,'Color',[.7 .7 .7]);
		else
			delete(handles.cfg.hcurdat(figIdx));
		end
	end
end
axes(handles.axes_xy);
handles.cfg.hcurdat(1) = plot(az,el,'k-','lineWidth',2);

axes(handles.axes_position);
handles.cfg.hcurdat(2) = plot(t,az,'b-','lineWidth',2);
hold on
handles.cfg.hcurdat(3) = plot(t,el,'r-','lineWidth',2);
title(handles.cfg.trial);

axes(handles.axes_velocity);
handles.cfg.hcurdat(4) = plot(t,vel,'g-','lineWidth',2);
hold on
handles.cfg.hcurdat(5) = plot(t,smv,'k-','lineWidth',2);
