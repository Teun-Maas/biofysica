close all
clearvars

d = '/Users/marcw/Dropbox/manuscripts/HRTF RECORDING DATA/Marielle HRTFs';

cd(d)

d = dir('*9996*.sphere');
fnames = {d.name};
nfiles = numel(fnames);

for ii = 1:nfiles
	dat = load(fnames{ii},'-mat');
	snd = dat.data.sndrec;
	mod = {dat.trialsingle.stim.modality}; 
	sel = strcmp(mod,'sound');
	az = dat.trialsingle.stim(sel).X;
	el = dat.trialsingle.stim(sel).Y;
	[az el]
	
	subplot(311)
	plot(snd(:,1)+ii/1000)
	hold on

	subplot(312)
	plot(snd(:,2)+ii/10000)
	hold on
	
	[f,m] = getpower(snd(:,2),48828.125);
	subplot(313)
	semilogx(f,m);
	hold on
	xlim([10 200])
	drawnow
end

