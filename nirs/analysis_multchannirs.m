function analysis_multchannirs(datFolder,dname,fname)

% basic NIRS analysis using fieldtrip functions
% requires the optodetemplates.xml in the added path
% addpath(genpath('D:\Matlab\PandA'))
% addpath(genpath('D:\Matlab\fieldtrip-20150421'))
% addpath('M:\Anja Roye\NIRS\Kennan experiment - 2 sides\LR-01-2015-06-16_passive')
% addpath('M:\Anja Roye\NIRS\Kennan experiment - 2 sides\analysis\')
% addpath('M:\Anja Roye\NIRS\Kennan experiment - 2 sides\analysis\fnirs.git\')

%% files
if nargin<3
	
	datFolder	= '/Volumes/mbaudit1/Marc van Wanrooij/NIRS/Kennan experiment - 2 sides';
	% 	fname		= 'LR-01-2015-06-16-0001_2.oxy3';
	
	% 	dname		= 'LR-AR-2015-06-01';
	% 	fname = 'LR-AR-2015-06-01-0002.oxy3';
	
	dname = 'LR-04-2015-06-17_active';
	fname = 'LR-04-2015-06-17-0001.oxy3';
	


% 		dname = 'LR-01-2015-06-16_passive';
% 		fname = 'LR-01-2015-06-16-0001_2.oxy3';
end
fname			= fullfile(datFolder,dname,fname);
rawfile			= fullfile(datFolder,dname,'data_raw.mat');
downfile		= fullfile(datFolder,dname,'data_down.mat');
transfile		= fullfile(datFolder,dname,'data_trans.mat');
triggerfile		= fullfile(datFolder,dname,'data_trigger.mat');
rcsfile			= fullfile(datFolder,dname,'data_rcs.mat');
stafile			= fullfile(datFolder,dname,'data_sta_epoched_rcs.mat');
devfile			= fullfile(datFolder,dname,'data_dev_epoched_rcs.mat');


avgstafile		= fullfile(datFolder,dname,'data_sta_avg_basel_rcs.mat');
avgdevfile		= fullfile(datFolder,dname,'data_dev_avg_basel_rcs.mat');

%% Flags
loadFlag		= [true true true true true true true];
% loadFlag		= [false false false false false false false];

%% Parameters
n				= 3150; % no. of samples to remove

%% Raw data
% read the raw data from file as one long continuous segment without any additional filtering
if ~loadFlag(1)
	cfg					= [];
	cfg.dataset			= fname;
	data_raw			= ft_preprocessing(cfg);
	save(rawfile, 'data_raw');
else
	load(rawfile);
end
% % Plot raw events
% % ISSUE: low-pass filtering of events
% figure(1)
% chansel  = [97 98];
% plot(data_raw.time{1}, data_raw.trial{1}(chansel(1), :)+4,'o-')
% hold on
% plot(data_raw.time{1}, data_raw.trial{1}(chansel(2), :),'o-')
% xlabel('time (s)')
% ylabel('channel amplitude (uV)')
% box off
% axis square
% legend(data_raw.label(chansel))


%% downsample
if ~loadFlag(2)
	
	cfg				= [];
	cfg.resamplefs	= 10;
	data_down		= ft_resampledata(cfg, data_raw);
	save(downfile, 'data_down');
else
	load(downfile);
end

% % Plot down-samples
% % ISSUE: what is event rate? 0.6 Hz, Seems close to heart-rate (simulation
% % of analysis needed)
% % ISSUE: small negative crosstalk between ADC channels
% figure(2)
% chansel  = [1 97 98];
% ax(1) = subplot(311);
% plot(data_down.time{1}, data_down.trial{1}(chansel(1), :),'Color',[1 0 1])
% hold on
% ax(2) = subplot(312);
% plot(data_down.time{1}, data_down.trial{1}(chansel(2), :))
% hold on
% ax(3) = subplot(313);
% plot(data_down.time{1}, data_down.trial{1}(chansel(3), :))
% xlabel('time (s)');
% linkaxes(ax,'x');
%
%
% x = data_down.trial{1}(chansel(2),:);
% figure;
% getpower(x',10,'display',true)


%%
%convert data, i.e. changes in optical density to concentration changes
if ~loadFlag(3)
	cfg = [];
	cfg.target	= {'O2Hb', 'HHb'};
	cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all?, all NIRS channels are selected then
	data_trans	= ft_transform_ODs(cfg, data_down);
	
	save(transfile, 'data_trans');
else
	load(transfile);
end

%% reference channel subtraction
% ISSUE: Does this work? Does this not introduce artefacts?
% Is heart rate removed?
if ~loadFlag(4)
	data_rcs	= nirs_rcs(data_trans,n);
	save(rcsfile, 'data_rcs');
	
else
	load(rcsfile)
end

% keyboard
% figure(3)
% chansel  = [1];
% plot(data_rcs.time{1}, data_rcs.trial{1}(chansel(1), :))
% xlabel('time (s)')

% ISSUE: heart rate removal via NAP

%%

%% get trigger indices based on raw data
% ISSUES: see above
% ISSUE: is number of triggers correct? 549 triggers in standard
if loadFlag(5)
	trig_chan		= [97 98];
	trgdur			= 100;
	fsdown			= data_down.fsample;
	data_trigger	= nirs_trigger(data_raw,trig_chan,trgdur,fsdown);
	save(triggerfile, 'data_trigger');
else
	load(triggerfile)
end

% pa_verline(data_rcs.time{1}(data_trigger.trigger{1}));
% drawnow

%%
on =  data_trigger.trigger{1};
off = on+0.1*10;
[nchan,N] = size(data_rcs.trial{1});
t = 1:N;
t = t/10;

X       = zeros(1, N);
for ii	= 1:length(on)
	X(on(ii):off(ii)) = 1;
end
Ystd = pa_nirs_hdrfunction([1 1 1],X);
close all
ax(1) = subplot(211);
plot(t,X,'k-','Color',[.7 .7 .7]);
hold on

plot(t,Ystd,'k-');
box off
plot(t,data_rcs.trial{1}(1,:))

on =  data_trigger.trigger{2};
off = on+0.1*10;
N = length(data_rcs.trial{1});

X       = zeros(1, N);
for ii	= 1:length(on)
	X(on(ii):off(ii)) = 1;
end
Ydev = pa_nirs_hdrfunction([1 1 1],X);
ax(2) = subplot(212);
plot(t,X,'k-','Color',[.7 .7 .7]);
hold on


plot(t,Ydev,'k-');
box off

linkaxes(ax,'x');


%% GLM
% ft_regressconfound
dat = data_rcs.trial{1};
x = [Ystd; Ydev];
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
if loadFlag(6)
	
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
else
	load(stafile);
	load(devfile);
end

%%


%%

%% Further topo and averaging

% keyboard

%% Averaging
if loadFlag(7)
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
else
	load(avgstafile);
	load(avgdevfile);
end


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
whos X Y Z	
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


%%
keyboard