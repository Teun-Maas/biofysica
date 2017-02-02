
%% Data folders
datFolder = '/Volumes/mbaudit1-1/Marc van Wanrooij/NIRS/Kennan experiment - 2 sides';
cd(datFolder);

d	= dir('LR*'); % search for experimenter LR
sel = [d.isdir];
d	= d(sel);

ndir = numel(d);
F = cell(ndir,1);
for ii = 2
	close all
	dname = d(ii).name
	
	cd(dname);
	f		= dir('*.oxy3');
	fname	= f.name
	
	
	nirs_analysis_glm_t(datFolder,dname,fname)
	cd ..
	
	F{ii} = fname;
	drawnow
end