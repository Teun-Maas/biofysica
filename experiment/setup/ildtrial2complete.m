function data = ildtrial2complete(dname)
% ILDTRIAL2COMPLETE
%
% Combine the ILD trial files to one file
%


%% Initialization
if nargin<1
	% 	fname = [];
	dir;
	d = what;
	dname = d.path;
end


cd(dname);
d = dir('*.ild');

%% find unique blocks
f		= char(d.name);
f		= f(:,1:end-12);
ublocks = unique(f,'rows');

%% Load date per block
nblocks = size(ublocks,1);
for blockIdx = 1:nblocks
	fname	= ublocks(blockIdx,:);
	d		= dir([fname '*.ild']);
	S = [];
	R = [];
	B = [];
	for fIdx = 1:numel(d)
		disp(['Loading ' d(fIdx).name])
		load(d(fIdx).name,'-mat');
		S	= [S stim]; %#ok<AGROW>
		R	= [R RT]; %#ok<AGROW>
		B	= [B button]; %#ok<AGROW>
		
	end
	
	data = [[S.ild]' [S.freq_left]' [S.freq_right]' B' R'];
	fname = fcheckext(fname,'ild');
	fname = fullfile(dname,fname);
	save(fname,'data','cfg');
end



