function handles = createdir(handles)
% [FNAME,DNAME] = GETFNAME(FPATH)
%
% Input dialog for data file name
%
% See also FEXISTDLG

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com

%% Initialization

cfg		= handles.cfg;

if ~exist(cfg.dname,'dir') % create directory
	mkdir(cfg.dname);
	mkdir([cfg.dname filesep 'trial']);
	str = ['Creating directory: ' cfg.dname];
	disp(str);
else
	str = {['   Directory ''' cfg.dname ''' already exists.'];'   Data file will be added to this folder.'};
	disp(char(str));
	mkdir([cfg.dname filesep 'trial']);
end

handles.cfg = cfg;
