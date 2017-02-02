function nirs_analysis_glm_t(datFolder,dname,fname)

% basic NIRS analysis using fieldtrip functions
% requires the optodetemplates.xml in the added path

%% files
if nargin<3
	datFolder	= '/Volumes/mbaudit1-1/Marc van Wanrooij/NIRS/Kennan experiment - 2 sides';
	dname = 'LR-01-2015-06-16_passive';
	fname = 'LR-01-2015-06-16-0001_2.oxy3';
end
fname			= fullfile(datFolder,dname,fname);
rawfile			= fullfile(datFolder,dname,'data_raw.mat');
downfile		= fullfile(datFolder,dname,'data_down.mat');
scifile			= fullfile(datFolder,dname,'data_sci.mat');
rejectfile		= fullfile(datFolder,dname,'data_reject.mat');
hemofile		= fullfile(datFolder,dname,'data_hemo.mat');
transfile		= fullfile(datFolder,dname,'data_trans.mat');
triggerfile		= fullfile(datFolder,dname,'data_trigger.mat');
rcsfile			= fullfile(datFolder,dname,'data_rcs.mat');
stafile			= fullfile(datFolder,dname,'data_sta_epoched_rcs.mat');
devfile			= fullfile(datFolder,dname,'data_dev_epoched_rcs.mat');
% avgstafile		= fullfile(datFolder,dname,'data_sta_avg_basel_rcs.mat');
% avgdevfile		= fullfile(datFolder,dname,'data_dev_avg_basel_rcs.mat');

%% Flags
close all;

%% Raw data
% read the raw data from file as one long continuous segment without any additional filtering
% data_raw = getraw(rawfile,fname);

%% Artifact rejection
% data_reject = getreject(rejectfile,data_raw); % to remove first samples, see also nirs_rmv
% data_reject = getreject(rejectfile); % to remove first samples, see also nirs_rmv

%% Hemoglobin
% data_hemo = gethemo(hemofile,data_reject);
% data_hemo = gethemo(hemofile);


%% downsample
% data_down = getdown(downfile,data_hemo);
% data_down = getdown(downfile);

%% Scalp coupling index
% data_sci = getsci(scifile,data_down);
% data_sci = getsci(scifile);

%% convert data, i.e. changes in optical density to concentration changes
% data_trans = gettrans(transfile,data_sci);
% data_trans = gettrans(transfile);

%% reference channel subtraction
% ISSUE: Does this work? Does this not introduce artefacts?
% Is heart rate removed?
% data_rcs = getrcs(rcsfile,data_trans);
data_rcs = getrcs(rcsfile);

% keyboard
%%



%% get trigger indices based on raw data
% ISSUES: see above
% ISSUE: is number of triggers correct? 549 triggers in standard
% 	load(triggerfile)

% trig_chan		= [97 98];
% trgdur			= 100;
% fsdown			= data_down.fsample;
% data_trigger	= nirs_trigger(data_raw,trig_chan,trgdur,fsdown);
% save(triggerfile, 'data_trigger');
% 
% %%
% on =  data_trigger.trigger{1};
% off = on+0.1*10;
% [nchan,N] = size(data_rcs.trial{1});
% t = 1:N;
% t = t/10;
% X       = zeros(1, N);
% for ii	= 1:length(on)
% 	X(on(ii):off(ii)) = 1;
% end
% Ystd = pa_nirs_hdrfunction([1 1 1],X);
% 
% figure(4)
% clf
% ax(1) = subplot(211);
% plot(t,X,'k-','Color',[.7 .7 .7]);
% hold on
% 
% plot(t,Ystd,'k-');
% box off
% plot(t,data_rcs.trial{1}(1,:))
% 
% on =  data_trigger.trigger{2};
% off = on+0.1*10;
% N = length(data_rcs.trial{1});
% 
% X       = zeros(1, N);
% for ii	= 1:length(on)
% 	X(on(ii):off(ii)) = 1;
% end
% Ydev = pa_nirs_hdrfunction([1 1 1],X);
% ax(2) = subplot(212);
% plot(t,X,'k-','Color',[.7 .7 .7]);
% hold on
% plot(t,Ydev,'k-');
% box off
% linkaxes(ax,'x');


%% GLM
% ft_regressconfound
dat = data_rcs.trial{1};
x = data_rcs.hemo{1};
nchan = size(dat,1);
T = NaN(nchan,2);
for ii = 1:nchan
	y			= dat(ii,:);
	stats		= regstats(y',x','linear','tstat');
	t			= stats.tstat.t(2:3);
	T(ii,:)		= t;
end

%% Segmenting continuous data into trials
%epoching
% ISSUE: High-pass/low-pass filter?

cfg						= [];
cfg.dataset				= fname;
cfg.trialdef.eventvalue = 1; % read conditions
cfg.trialdef.prestim    = 5; % in seconds
cfg.trialdef.poststim   = 15; % in seconds

sta						= data_trigger.trigger{1}; % trigger indices for onset of the standard
dev						= data_trigger.trigger{2}; % trigger indices for onset of the deviant

trl						= NaN(length(sta)-2,4);
for stdidx = 1:length(sta)-2
	begsample			= sta(stdidx) - cfg.trialdef.prestim*data_down.fsample;
	endsample			= sta(stdidx) + cfg.trialdef.poststim*data_down.fsample - 1;
	offset			= -cfg.trialdef.prestim*data_down.fsample;
	trigger			= 1; % remember the trigger (=condition) for each trial
	trl(stdidx, :)		= [round([begsample endsample offset])  trigger];
end

cfg.trl			= trl;
cfg.lpfilter	= 'yes';                              % apply lowpass filter
cfg.lpfreq		= 10;                                 % lowpass
%cfg.hpfilter   = 'yes';                              % apply highpass filter
%cfg.hpfreq     = 0.02;                                 %
data_sta_epoched = ft_redefinetrial(cfg,data_rcs);

%% deviant epoching
cfg						= [];
%cfg.dataset  = 'LR-AR-2015-04-21-0001.oxy3';
cfg.dataset				= fname;
cfg.trialdef.eventvalue = 2; % read conditions
cfg.trialdef.prestim    = 5; % in seconds
cfg.trialdef.poststim   = 15; % in seconds

trl						= NaN(length(dev),4);
for stdidx=1:length(dev)
	begsample     = dev(stdidx) - cfg.trialdef.prestim*data_down.fsample;
	endsample     = dev(stdidx) + cfg.trialdef.poststim*data_down.fsample - 1;
	offset        = -cfg.trialdef.prestim*data_down.fsample;
	trigger       = 2; % remember the trigger (=condition) for each trial
	trl(stdidx, :)	= [round([begsample endsample offset])  trigger];
end
cfg.trl = trl;

cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 10;                                 % lowpass
%cfg.hpfilter   = 'yes';                              % apply highpass filter
%cfg.hpfreq     = 0.02;                                 %
data_dev_epoched = ft_redefinetrial(cfg,data_rcs);


save(stafile, 'data_sta_epoched');
save(devfile, 'data_dev_epoched');


%%


%%

%% Further topo and averaging

% keyboard

%% Averaging
% else
% 	load(avgstafile);
% 	load(avgdevfile);
% end

cfg = [];

%cfg.channel = 'gui'; % opens a gui to select channels
data_sta_avg_rcs = ft_timelockanalysis(cfg, data_sta_epoched);
% save([DataFolder 'data_sta_avg_rcs.mat'], 'data_sta_avg_rcs');
data_dev_avg_rcs = ft_timelockanalysis(cfg, data_dev_epoched);
% save([DataFolder 'data_dev_avg_rcs.mat'], 'data_dev_avg_rcs');


%to perform baseline normalization
cfg.baseline = [-5 0];
data_sta_avg_basel_rcs = ft_timelockbaseline(cfg, data_sta_avg_rcs);
% save([DataFolder 'data_sta_avg_basel_rcs.mat'], 'data_sta_avg_basel_rcs');
data_dev_avg_basel_rcs = ft_timelockbaseline(cfg, data_dev_avg_rcs);
% save([DataFolder 'data_dev_avg_basel_rcs.mat'], 'data_dev_avg_basel_rcs');
%save timelock



%%
%% plotting
% ISSUE: How is layout defined?
% cfg = [];
% cfg.image = 'optotemp.png';
% lay = ft_prepare_layout(cfg);
%
% cfg.layout = lay;
% ft_layoutplot(cfg, data_trigger)
load('lay.mat')

figure(100)
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = lay;
cfg.interactive = 'yes';
cfg.channel = '* [O2Hb]';
cfg.ylim = [-0.2 0.2];
%cfg.graphcolor = 'r';
ft_multiplotER(cfg, data_sta_avg_basel_rcs,data_dev_avg_basel_rcs)

figure(101)
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = lay;
cfg.interactive = 'yes';
cfg.channel = '* [HHb]';
cfg.ylim = [-0.2 0.2];
%cfg.graphcolor = 'b';
ft_multiplotER(cfg,  data_sta_avg_basel_rcs,data_dev_avg_basel_rcs);

%%
% keyboard

%% topolot
% ISSUE: topoplot - layout should be recreate, outline - sideview
% for top-view: cf standard lay out files?
% close all
cfg				= [];
cfg.xlim		= [5 10];
cfg.channel = '* [O2Hb]';

figure
cfg.layout		= 'kennanhelmet.lay';
ft_layoutplot(cfg)

% figure(2)
% cfg.layout		= 'EEG1020.lay';
% ft_layoutplot(cfg)
%
%
% figure(3)
% cfg.layout		= 'EEG1010.lay';
% ft_layoutplot(cfg)

%%
% M = [[1:66]' lay.pos lay.width lay.height];
% xlswrite('kennannirs.xlsx',lay.label)
% dlmwrite('kennannirs.lay',M,'delimiter','\t')
%
cfg.colorbar	= 'yes';
cfg.channel = '* [O2Hb]';

data = data_dev_avg_basel_rcs;
data.avg = repmat(T(:,2),1,200);

ft_topoplotER(cfg,data);

% ax = [-0.6250    0.6250   -0.5076    0.5750];
% h = patch([.2*ax(1) .2*ax(1) .2*ax(2) .2*ax(2)],[ax(3) ax(4) ax(4) ax(3)],[.95 .95 .95]);
% get(h)
% set(h,'EdgeColor','none');

%% Interpolation
figure(666)
stdidx = [1 1 2 2];
side = [1 2 1 2];
for ii = 1:4
	switch side(ii)
		case 1
			[Xq,Yq] = meshgrid(-30:-10,-10:10);
		case 2
			[Xq,Yq] = meshgrid(10:30,-10:10);
			
	end
	[X,Y,~] = nirs_optodelayout(data_trans,'disp',false,'lay',false);
	% 	X		= X;
	% 	Y		= Y;
	Z		= T(1:2:end,stdidx(ii));
	switch side(ii)
		case 1
			sel		= X<0;
		case 2
			sel = X>0;
	end
	X		= X(sel);
	Y		= Y(sel);
	Z		= Z(sel);
	
	
	
	% 	F		= scatteredInterpolant(X,Y,Z,'natural','none');
	% 	Vq		= F(Xq,Yq);
	
	Vq = griddata(X,Y,Z,Xq,Yq,'cubic');
	subplot(2,2,ii)
	contourf(Xq,Yq,Vq,-20:5:20)
	colorbar;
	axis square
	switch stdidx(ii)
		case 1
			title('Standard');
		case 2
			title('Deviant');
			
	end
	shading flat
end

C = [];
for ii = 1:4
	subplot(2,2,ii)
	cax = caxis;
	C = [C;cax];
end
mx = max(abs(C(:)));
for ii = 1:4
	subplot(2,2,ii)
	caxis([-mx mx]);
end



function data_raw = getraw(rawfile,fname)

tic
if exist(rawfile,'file')
	disp('Loading raw');
	load(rawfile);
	% 		if ~isfield(data_raw,'nrm'); %#ok<NODEF>
	% 			data_raw		= nirs_rmvsamples(data_raw); % remove first samples
	% 			save(rawfile, 'data_raw');
	% 		end
	
else
	cfg					= [];
	cfg.dataset			= fname;
	cdf.demean			= 'yes';
	data_raw			= ft_preprocessing(cfg);
	% 	data_raw		= nirs_rmvsamples(data_raw);
	save(rawfile, 'data_raw');
end
toc


% Plot raw events
% ISSUE: low-pass filtering of events
dat = data_raw.trial{1};
nchan = size(dat,1);
sb = ceil(sqrt(nchan));

for ii = 1:nchan
	figure(1)
	subplot(sb,sb,ii);
	plot(dat(ii,:))
	hold on
	
	box off;
	axis off;
	% 		ylim([0 6.5]);
	xlim([0 length(dat(ii,:))]);
	str = [data_raw.label{ii}];
	title(str)
end


figure(2)
clf
chansel		= [1 97 98];
ax(1)		= subplot(343);
plot(data_raw.time{1}, data_raw.trial{1}(1:96, :))
hold on
title('Channel 1');

ax(2)		= subplot(347);
plot(data_raw.time{1}, data_raw.trial{1}(chansel(2), :))
hold on
xlabel('time (s)');
title('Standard event');

ax(2)		= subplot(3,4,11);
plot(data_raw.time{1}, data_raw.trial{1}(chansel(3), :))
hold on
xlabel('time (s)');
title('Standard event');

x = data_raw.trial{1}(chansel(1),:);
subplot(344)
getpower(x',10,'display',true);
set(gca,'XTick',[0.1 1 10],'XTickLabel',[0.1 1 10]);

x = data_raw.trial{1}(chansel(2),:);
subplot(348)
getpower(x',10,'display',true);
set(gca,'XTick',[0.1 1 10],'XTickLabel',[0.1 1 10]);


x = data_raw.trial{1}(chansel(3),:);
subplot(3,4,12)
getpower(x',10,'display',true);
set(gca,'XTick',[0.1 1 10],'XTickLabel',[0.1 1 10]);

function data_down = getdown(downfile,data)
tic
if exist(downfile,'file')
	load(downfile);
else
	
	cfg				= [];
	cfg.resamplefs	= 10;
	cfg.demean		= 'yes';
	% downsample trial data
	data_down		= ft_resampledata(cfg, data);
	
	% downsample hemodynamic prediction
	data.trial{1} = data.hemo{1};
	data.label = data.hemolabel;

	data_hemo		= ft_resampledata(cfg, data);
	data_down.hemo{1} = data_hemo.trial{1};
	
	% save
	save(downfile, 'data_down');
end
toc

%% Plot down-samples
% ISSUE: what is event rate? 0.6 Hz, Seems close to heart-rate (simulation
% of analysis needed)
% ISSUE: small negative crosstalk between ADC channels
tic

%%
figure(2)
chansel		= [1 97 98];
ax(1)		= subplot(341);
plot(data_down.time{1}, data_down.trial{1}(1:96, :))
hold on
title('Channel 1');

ax(2)		= subplot(345);
plot(data_down.time{1}, data_down.trial{1}(chansel(2), :))
hold on
xlabel('time (s)');
title('Standard event');

ax(2)		= subplot(349);
plot(data_down.time{1}, data_down.trial{1}(chansel(3), :))
hold on
xlabel('time (s)');
title('Standard event');

x = data_down.trial{1}(chansel(1),:);
subplot(342)
getpower(x',10,'display',true);
set(gca,'XTick',[0.1 1 10],'XTickLabel',[0.1 1 10]);

x = data_down.trial{1}(chansel(2),:);
subplot(346)
getpower(x',10,'display',true);
set(gca,'XTick',[0.1 1 10],'XTickLabel',[0.1 1 10]);


x = data_down.trial{1}(chansel(3),:);
subplot(3,4,10)
getpower(x',10,'display',true);
set(gca,'XTick',[0.1 1 10],'XTickLabel',[0.1 1 10]);
toc

function data_reject = getreject(rejectfile,data_raw)

if exist(rejectfile,'file')
	load(rejectfile);
else
	
	load('lay.mat')
	cfg = [];
	% cfg.showlabels = 'yes';
	cfg.layout = lay;
	cfg.method   = 'channel';
	% 	cfg.channel = '* [*]';
	
	% data_reject = ft_rejectvisual(cfg,data_raw);
	cfg						= ft_databrowser(cfg,data_raw);
	cfg.artfctdef.reject	= 'partial';
	data_reject				= ft_rejectartifact(cfg,data_raw);
	save(rejectfile,'data_reject');
end

function data_hemo = gethemo(hemofile,data)

if exist(hemofile)
	load(hemofile);
else
	chansel = [97 98]; % ADC, 1 = std, 2 = dev
	
	X = data.trial{1}(chansel(1),:);
	Ystd = pa_nirs_hdrfunction([1 1 1],X,'Fs',data.fsample);
	
	X = data.trial{1}(chansel(2),:);
	Ydev = pa_nirs_hdrfunction([1 1 1],X,'Fs',data.fsample);
	
	dat = data.trial{1};
	n	= size(dat,1);
	
	data.hemo{1} = [Ystd; Ydev];
	data.hemolabel{1} = 'Standard Hemo';
	data.hemolabel{2} = 'Deviant Hemo';
	
	data_hemo = data;
	save(hemofile,'data_hemo');
end

function data_sci = getsci(scifile,data_down)
if exist(scifile,'file')
	load(scifile);
else
data_sci = nirs_sci(data_down,'disp',true);
save(scifile, 'data_sci');
end

function data_trans = gettrans(transfile,data_sci)

tic
if exist(transfile,'file');
	load(transfile);
else
	cfg = [];
	cfg.target	= {'O2Hb', 'HHb'};
	cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all?, all NIRS channels are selected then
	data_trans	= ft_transform_ODs(cfg, data_sci);
	
	%% select bad channels
	sci = data_trans.sci;
	sel = sci>0.75;

	data_trans.label			= data_trans.label(sel);
	data_trans.trial{1}			= data_trans.trial{1}(sel,:);
	data_trans.opto.chanpos		= data_trans.opto.chanpos(sel);
	data_trans.opto.chantype	= data_trans.opto.chantype(sel);
	data_trans.opto.label		= data_trans.opto.label(sel);
	data_trans.opto.DPF			= data_trans.opto.DPF(sel);
	
	%% Determine optodedistance
	data_trans					= nirs_optodedistance(data_trans);
	save(transfile, 'data_trans');
	
end
toc

function data_rcs = getrcs(rcsfile,data_trans);
figure(6)
if exist(rcsfile,'file');
	load(rcsfile)
	disp('Load rcs file')
else
	data_rcs	= nirs_rcs(data_trans);
	save(rcsfile, 'data_rcs');
end
