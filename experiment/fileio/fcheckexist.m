function fname = pa_fcheckexist(fname,fspec)
% FNAME = PA_FCHECKEXIST(FNAME)
%
% Check whether file FNAME exists. If not, explorer window will open.
% Output filename FNAME will be from the chosen file.
%
% FNAME = PA_FCHECKEXIST(FNAME,FSPECI) allows to add a file specifier (by
% default: '*.*');
%
% See also PA_FCHECKEXT

% (c) 2011 dr. Marc M. van Wanrooij

%% Initialization
if nargin<2
    fspec               = '*.*'; % check all files
end
if nargin<1
    fname               = ''; % without input, PA_FCHECKEXIST is useful for opening an explorer window
end
if isempty(fname)
    str                 = 'Choose a file';
else
    [pathstr,name,ext]  = fileparts(fname); %#ok<ASGLU>
    str                 = ['[' name ext '] not found, Choose a file'];
end

%% Open explorer window
if ~exist(fname,'file') % When file does not exist, 
    [fname,pname]       = uigetfile(fspec,str); % open explorer window
    fname               = [pname fname];
    if ~ischar(fname) % cancel has been pressed
        fname           = '';
		disp('No filename has been chosen. This will cause errors');
        return
    end
end
