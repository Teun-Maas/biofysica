function cleanspherefiles(dname,EI)
% Remove graphic handles from sphere-files

if nargin<1
	dname = cd;
end
cd(dname)

if nargin<2
	d = dir('*-*-*-*');
	fname = d(1).name;
	EI = fname(1:2);
end

d = dir([EI '-*']);
ndir = numel(d);
for dIdx = 1:ndir
	dname = d(dIdx).name;
	cd(dname)
	cd('trial')
	f = dir([EI '-*.sphere']);
	nfiles = numel(f);
	for fIdx = 1:nfiles
		fname = f(fIdx).name;
		load(fname,'-mat')
		if isfield(cfg,'hcurtar')
			disp(fname)
			cfg = rmfield(cfg,'hcurtar');
		end
		if isfield(cfg,'hcurdat')
			disp(fname)
			cfg = rmfield(cfg,'hcurdat');
		end
		save(fname,'data','trialsingle','cfg','-mat');
		% 	whos
	end
	
	cd ..
	cd ..
end