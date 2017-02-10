function handles = gethandles(handles)
cfg					= handles.cfg;
expIdx				= get(handles.popupmenu_expinitials,'Value');
expInitials			= get(handles.popupmenu_expinitials,'String');
cfg.expInitials		= expInitials{expIdx};
cfg.subjectid		= get(handles.edit_subjectid,'String');
cfg.datestring		= get(handles.edit_date,'String');
cfg.block			= get(handles.edit_block,'String');
cfg.fname			= [cfg.expInitials '-' sprintf('%04u',str2double(cfg.subjectid)) '-' cfg.datestring '-' sprintf('%04u',str2double(cfg.block)) '.sphere']; % file name
cfg.dname			= [cfg.fpath cfg.expInitials filesep cfg.expInitials '-' sprintf('%04u',str2double(cfg.subjectid)) '-' cfg.datestring]; % directory name
cfg.expdir			= [cfg.fpath cfg.expInitials filesep 'EXP' filesep]; % exp directory name
cfg.snddir			= [cfg.fpath cfg.expInitials filesep 'SND' filesep]; % wav directory name
expIdx				= get(handles.popupmenu_expinitials,'Value');

% get exp and cfg files
str					= [cfg.expdir filesep '*.exp'];
d					= dir(str); % default exp folder
cfg.expfiles		= {d.name};
if isempty(cfg.expfiles)
	
end
set(handles.popupmenu_exp,'String',cfg.expfiles)
expfileIdx			= get(handles.popupmenu_exp,'Value');
cfg.expfname		= cfg.expfiles{expfileIdx};

d					= dir([cfg.expdir filesep '*.cfg']);
cfg.cfgfiles		= {d.name};
set(handles.popupmenu_cfg,'String',cfg.cfgfiles);
cfgfileIdx			= get(handles.popupmenu_cfg,'Value');
cfg.cfgfname		= cfg.cfgfiles{cfgfileIdx};
handles.cfg			= cfg;
