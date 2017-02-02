function fname = pa_fcheckext(fname,fext)
% FNAME = PA_FCHECKEXT(FNAME,FEXT)
%
% Check whether extension of file FNAME corresponds to FEXT
% If not, the extension will be replaced or added.
%
% See also PA_FCHECKEXIST

% (c) 2011 Marc van Wanrooij

%% Check dot
if ~strcmp(fext(1),'.')
    fext = ['.' fext];
end;

[pathstr,name,ext] = fileparts(fname);

%% Check extension
if ~strcmp(ext,fext)
    ext     = fext;
end

%% Reconstruct filename
fname       = fullfile(pathstr,[name ext]);
