function dirname = pa_datadir(expname)
% PA_DATADIR
% Jump to the data-directory:
% C:\DATA
%
% It is useful for groups of people working on the same data-sets to use the
% same directory on their computer for data storage. This function might
% help even more.
%
% PA_DATADIR('EXP')
% Jump to the data-directory:
% C:\DATA\EXP
%
% Fileseparator is determined for OS.
%
% see also CD, FILESEP

% (c) 2011-05-06 Marc van Wanrooij

F		= filesep;
dirname = ['C:' F 'DATA' F];
if ~exist(dirname,'dir')
% 	dirname = ['E:' F 'DATA' F];
	dirname = ['/Users/marcw/DATA'];
end
if nargin>0
	dirname = [dirname expname];
end

if exist(dirname,'dir') % check if directory exists
    cd(dirname)
else % or open a GUI
    dirname = uigetdir;
    if ~ischar(dirname)
        return
    end
    cd(dirname)
end
