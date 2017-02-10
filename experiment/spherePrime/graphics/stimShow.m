function handles = stimShow(handles)

axes(handles.axes_xy); %#ok<*NASGU>
if isfield(handles.cfg,'hcurtar')
	n = numel(handles.cfg.hcurtar);
	for ii = 1:n
		if isfield(get(handles.cfg.hcurtar(ii)),'MarkerSize')
			set(handles.cfg.hcurtar(ii),'MarkerSize',5,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[.7 .7 .7]);
		end
	end
end
for stmIdx = 1:numel(handles.stim)
	switch upper(handles.stim(stmIdx).modality)
		case {'LED','SKY'}
			handles.cfg.hcurtar(stmIdx) = plot(handles.stim(stmIdx).azimuth,handles.stim(stmIdx).elevation,'r*','MarkerSize',8);
		case {'SOUND','SND','SND1','SND2','SND3','SND4'}
			handles.cfg.hcurtar(stmIdx) = plot(handles.stim(stmIdx).azimuth,handles.stim(stmIdx).elevation,'bo','MarkerSize',8);
	end
end
