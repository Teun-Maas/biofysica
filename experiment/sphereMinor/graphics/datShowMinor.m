function handles = datShowMinor(handles)


%% load data
Hchan		= 1; Vchan		= 2; Fchan		= 3; % set channels in correct order
S			= load(handles.cfg.calfile,'-mat'); % Check: correct for sphereMinor?
dat			=  handles.data(1).raw;
H			= dat(:,Hchan); H = H(:);
V			= dat(:,Vchan); V = V(:);
F			= dat(:,Fchan); F = F(:);

cfactor		= 1.7e4; % conversion factor Volt --> deg
H			= H*cfactor;
V			= V*cfactor;
F			= F*cfactor;

DAT			= [H V F]';
%% calibrate
[az,el]		= pa_calib(DAT,S); % calibrate: from Volts --> Deg (correct for sphereMinor?)
% az          = filtfilt(handles.cfg.lpFilt,az); % filter
% el          = filtfilt(handles.cfg.lpFilt,el); % filter


t			= (1:numel(H))/handles.cfg.RZ6Fs; % time (s)
[vel,smv]	= getvel(az',el',handles.cfg.RZ6Fs,0.01); % get velocities

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
handles.cfg.hcurdat(1) = plot(H,V,'k-','lineWidth',2);

axes(handles.axes_position);
handles.cfg.hcurdat(2) = plot(t,H,'b-','lineWidth',2);
hold on
handles.cfg.hcurdat(3) = plot(t,F,'r-','lineWidth',2);
title(handles.cfg.trial);

axes(handles.axes_velocity);
handles.cfg.hcurdat(4) = plot(t,vel,'g-','lineWidth',2);
hold on
handles.cfg.hcurdat(5) = plot(t,smv,'k-','lineWidth',2);


%% provisional: plot target/response curves with target to final position
% there's probably a smarter way to do this, but this is easiest.
selsnd		= strcmpi({handles.stim.modality},'sound');
if any(selsnd)
	azTarget = handles.stim(selsnd).azimuth;
	elTarget = handles.stim(selsnd).elevation;

	axes(handles.axes_regress_azimuth);
	handles.cfg.hcurdat(6) = plot(azTarget,H(end),'b.', 'MarkerSize',12);
	axes(handles.axes_regress_elevation);
	handles.cfg.hcurdat(7) = plot(elTarget,V(end),'b.', 'MarkerSize',12);
	
end