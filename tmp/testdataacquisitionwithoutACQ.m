clearvars
close all

cd('/Users/marcw/DATA/Katharina Vogt/KV-0020-16-05-12/trial');

d = dir('*-0001-0*.sphere');

nfiles = numel(d);
n=0;
col = hot(6);
% for ii = 1:nfiles
	for ii = 1:9
	fname = d(ii).name;
	load(fname,'-mat');
	D = data.raw;
	
	disp('----------')
	t = trialsingle.stim;
	mod = {t.modality};
	if any(strcmp('data acquisition',mod))
		figure
		n = 0;
	end
	n = n+1;
	subplot(121)
	plot(D(:,1)+n*0.1,'Color',col(n,:));
	hold on
	ylim([-10 10]);
	title(n)
	subplot(122)
	plot(D(:,2)+n*0.1,'Color',col(n,:));
	hold on
	ylim([-10 10]);
	
	title(mod)
end