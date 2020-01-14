clearvars
close all


d = dir('._*')
fnames = {d.name}
nfiles = numel(fnames);

for ii = 1:nfiles
	fname = fnames{ii}
	delete(fname)
end