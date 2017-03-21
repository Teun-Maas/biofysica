function cleanspherefiles(exp)
% CLEANSPHEREFILES(DIR,EXP)
%
% Cleaning the .sphere files by removing graphic handles in the current
% directory of experimenter EXP
%

d = dir([exp '-*']);
ndir = numel(d);
for dIdx = 1:ndir
	dname = d(dIdx).name;
	cd(dname)
	cd('trial')
	f	= dir('SS-*.sphere');
	nfiles = numel(f);
	for fIdx = 1:nfiles
		fname = f(fIdx).name;
		load(fname,'-mat')
		if isfield(cfg,'hcurtar')
			cfg = rmfield(cfg,'hcurtar');
		end
		if isfield(cfg,'hcurdat')
			cfg = rmfield(cfg,'hcurdat');
		end
		save(fname,'data','cfg','-mat');
	end
	cd ..
	cd ..
end