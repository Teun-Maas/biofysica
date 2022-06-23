function S = loadmatindir(dname,str)
% S = LOADMATINDIR
%
% loads all mat files in current directory, and saves data in structure S
%
% S = LOADMATINDIR(DNAME)
%
% loads all mat files in directory DNAME
%
% S = LOADMATINDIR(DNAME,STR)
%
% loads all mat files with STR at the start of the filename in directory DNAME
%
% See also LOAD, DIR, CD

%% Initialization
if nargin<1
	dname = cd;
end
if nargin<2
	str = '';
end
cd(dname);


%% Determine mat files
d		= dir([ str '*.mat']);
fnames	= {d.name};
nfiles	= numel(fnames);
for ii	= 1:nfiles
	fname = fnames{ii};
	s = load(fname);
	S(ii).s = s;
end

%% End




