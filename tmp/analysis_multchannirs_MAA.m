% function tmp(fname,d)
% Use tmp.m as temporary script/batch function
function analysis_multchannirs_MAA

close all;
clearvars;

%% Checks
checkTrigger	= false;
checkClean		= false;

%% Initialization
d		= '/Users/marcw/DATA/Roos/NIRS/RC-CI-AV-05-2016-04-06'; % data directory
fname	= 'RC-CI-AV-05-2016-04-06-0001.oxy3'; % file name

d	=	'/Users/marcw/DATA/Peter Bremen/NIRS/Data/PB/Session2/';	%-- 6000 ms stimuli --%
fname	=	'MAA_Pilot_PB_MW_04_05_2016.oxy3';

% d	=	'/Users/marcw/DATA/Peter Bremen/NIRS/Data/PB/Session1/';	%-- 100 ms stimuli --%
% fname	=	'MAA_Pilot_MW_PB_2016_04_29.oxy3';
% 
d	=	'/Users/marcw/DATA/Peter Bremen/NIRS/Data/MW/Session1/';	%-- 6000 ms stimuli --%
fname	=	'MAA_Pilot_PB_MW_05_09_2016.oxy3';




%%
% Go to data directory
cd(d); % goto

header = read_artinis_oxy3(fname);
data   = read_artinis_oxy3(fname, header);

%% Data
[data_raw,StmMtx]		= loaddata(d,fname);
data_trigger	= gettime(data_raw,checkTrigger); % Timing / synchronisation
data_clean		= cleandata(data_trigger,checkClean); % remove on and offset
data_down		= resampledata(data_clean); % Resample
data_sci		= nirs_sci(data_down,'disp',false); % Scalp coupling index
data_nophys		= nirs_removephysnoise(data_sci); % remove heartbeat
data_trans		= transformdata(data_nophys); % transformation
% data_tmp		= nirs_rcs(data_trans); % Reference channel subtraction
data_rcs		= getdeep(data_trans); % only deep channels
side = 2; % 1 or 2: left vs right ?
data_epoched	= epochdata(data_rcs,StmMtx,side); % create epochs


%% Visualization

trial	= data_trans.trial{1};
t		= data_trans.time{1}/60; % min

figure;
for trlnr = 1:2:size(trial,1)
	plot(t,trial(trlnr,:));
	hold on
	str = data_trans.label{trlnr};
	% pause(0.5);
	% drawnow
end
title(str);
xlabel('Time (min)');

%%
trial	= data_rcs.trial{1};
t		= data_rcs.time{1}/60; % min
sync	= data_rcs.sync{1};
figure
nchan	= size(trial,1);
n		= ceil(sqrt(nchan));
for trlnr = 1:size(trial,1)
	ax(trlnr) = subplot(n,n,trlnr);
	plot(t,zscore(trial(trlnr,:)));
	hold on
	plot(t,zscore(sync))
	% 	plot(t,trial(trlnr,:));
	str = data_rcs.label{trlnr};
	title(str)
	ylim([-1 1]);
	% pause(0.5);
	% drawnow
end
linkaxes(ax,'x');
xlabel('Time (min)');


%%
cfg				= [];
cfg.baseline	= [-5 0];
data_avg		= ft_timelockanalysis(cfg, data_epoched);

%%
load('lay.mat')

figure(100)
cfg				= [];
cfg.showlabels	= 'yes';
cfg.layout		= lay;
cfg.interactive = 'yes';
cfg.channel		= '* [O2Hb]';
cfg.ylim		= [-0.2 0.2];
ft_multiplotER(cfg, data_avg)


%%

%% GLM
% ft_regressconfound
% dat = data_rcs.trial{1};
% x = data_rcs.hemo{1};
% nchan = size(dat,1);
% T = NaN(nchan,2);
% for ii = 1:nchan
% 	y			= dat(ii,:);
% 	stats		= regstats(y',x','linear','tstat');
% 	t			= stats.tstat.t(2:3);
% 	T(ii,:)		= t;
% end

%%
% figure
% cfg.layout		= 'kennanhelmet.lay';
% ft_layoutplot(cfg)

% %%
% close all
% plot(data_avg.time,data_avg.avg')
% data_avg.time

% %%
% cfg = [];
% cfg.colorbar	= 'yes';
% cfg.channel = '* [O2Hb]';
% cfg.layout		= lay;
%
% cfg.layout		= 'kennanhelmet.lay';
%
% d = data_avg.avg;
%
% sel = isnan(d);
% sel = sum(sel)>1;
% whos d
% d = d(:,sel);
% data = data_avg;
% data.avg = d;
% whos d
% ft_topoplotER(cfg,data);
keyboard
%%
T = nanmax(data_avg.avg,[],2);

[X,Y,~] = nirs_optodelayout(data_trans,'disp',[false false],'lay',false);
Z = T(1:2:end);
f = figure;
side = [1 2];
for ii = 1:2
	switch side(ii)
		case 1
			[Xq,Yq] = meshgrid(-30:-10,-10:10);
			sel		= X<0;
			
		case 2
			[Xq,Yq] = meshgrid(10:30,-10:10);
			sel = X>0;
			
			
	end
	
	whos sel X Y Z Xq Yq
	Vq = griddata(X(sel),Y(sel),Z(sel),Xq,Yq,'cubic');
	% 	subplot(2,2,ii)
	contourf(Xq,Yq,Vq)
	hold on
	% 	colorbar;
	% 	axis square
	% 	axis off
	cax = caxis;
	c = max(abs(cax));
	caxis([-c c])
end
[X,Y,~] = nirs_optodelayout(data_trans,'disp',[false true],'lay',false,'fig',f);
axis equal
%%
keyboard



%%

function plottrig(data_trigger)
figure
trial	= data_trigger.trial{1};
t		= data_trigger.time{1}/60; % min

sync		= 	trial(97:100,:);
% 1 = movie
% 2 = movie
% 3 = user interface
% 4 = nothing (screeen refresh)
% remove glitches from screen refresh
b			= regstats(sync(3,:),sync(4,:),'linear','r');
sync(3,:)	= b.r;
sync(3,:)	= sync(3,:)-sync(3,1); % offset

% graphic check
figure(2)
ax1(1) = subplot(211);
plot(sync(3,:));
ax1(2) = subplot(212);
plot(trial(99,:));
linkaxes(ax1,'x');

% then 'digitize' response
threshold = 0.15; % arbitrary, should be low enough to include onset, but large enough not to include glitches
sync(3,:) = sync(3,:)>threshold;

% graphic check
subplot(211);
hold on
plot(sync(3,:));

% Detect on- and offset
% use differences
whos sync

d		= [0 diff(sync(3,:))];
onset	= find(d>0);
offset	= find(d<0);

% How many?
nmovie = numel(onset); % should be 120
nstim = 120;
if nmovie~=nstim
	error('NMovieDoesNotMatchNStim');
else
	disp([num2str(nstim) ' movies were presented']);
end
% graphic check
figure(3)
ax1(1) = subplot(221);
plot(sync(3,:));
ax1(2) = subplot(222);
plot(d);
subplot(223);
hist(diff(t(onset)*60),0:5:500);
linkaxes(ax1,'x');

%
sync = sync(2:3,:);

figure(4)
plot(sync(1,:));

function Tevent = geteventtrig(Trg,Time,Fs,Show)

if( nargin < 4 )
	Show	=	0;
end

Trg		=	Trg ./ max(Trg);

[K,X]	=	genkernel(0,.01,Fs);

cTrg	=	conv(Trg,K,'same');
cTrg	=	cTrg ./ max(cTrg);

[~,Loc]	=	findpeaks(abs(cTrg),'MINPEAKHEIGHT',0.5);
Pks		=	cTrg(Loc);
Npks	=	length(Pks);

whos Time Loc
DetPtime=	Time(Loc);

Tevent	=	reshape(DetPtime,2,Npks/2)';

if( Show )
	xx		=	[min(Time) max(Time)];
	yy		=	[min(cTrg)*1.1 max(cTrg)*1.1];
	
	subplot(1,2,1)
	plot(X,K,'k-')
	xlim([min(X) max(X)])
	ylim([min(K)*1.1 max(K)*1.1])
	set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
	xlabel('time [s]')
	ylabel('amplitude [a.u.]')
	title('Gaussian kernel')
	axis('square')
	
	subplot(1,2,2)
	plot(Time,Trg,'k-')
	hold on
	plot(Time,cTrg,'r-')
	for k=1:Npks
		plot(DetPtime(k),Pks(k),'ko','LineWidth',1,'MarkerSize',8)
	end
	xlim(xx)
	ylim(yy)
	set(gca,'FontName','Arial','FontWeight','bold','FontSize',12)
	xlabel('time [s]')
	ylabel('amplitude [a.u.]')
	title('pulse with edge detection')
	axis('square')
end

function [k,X] = genkernel(mu,si,Fs)

X	=	-10:1/Fs:10;
G	=	exp( -( (X-mu).^2 ./ (2*si.^2) ) );
dG	=	gradient(G);

k	=	dG ./ max(dG);

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

function [data_raw,StmMtx] = loaddata(d,fname)

load([d fname(1:end-5) '.mat']);
StmMtx	=	Dat;
clear Dat

% read the raw data from file as one long continuous segment without any additional filtering
cfg					= [];
cfg.dataset			= fname;
data_raw			= ft_preprocessing(cfg);


function data_trigger = gettime(data_raw,checkTrigger)
%-- Calculate the time axis --%
Time					= data_raw.time{1};
sel						= strcmp(data_raw.label,'ADC005');
Trig					= data_raw.trial{1}(sel,:);
%-- Get event timestamps --%
Tevent					=	geteventtrig(Trig,Time,data_raw.fsample);

data_trigger				= data_raw;
data_trigger.trigger		= {Tevent*data_raw.fsample};
data_trigger.triggertime	= {Tevent};
data_trigger.sync = {Trig};


sync			= data_trigger.sync{1};
triggertime		= data_trigger.triggertime{1};
time			= data_trigger.time{1};

figure
plot(time,sync)
verline(triggertime)

if checkTrigger
	plottrig(data_trigger)
end

function data_clean = cleandata(data_trigger,checkClean)

a					= [data_trigger.trigger{1}];
t					= data_trigger.time{1}; % s
selt				= t>a(1)/data_trigger.fsample-10 & t<a(end)/data_trigger.fsample+10; % within 10 s of on- or offset
data_clean			= data_trigger;
b					= data_clean.time{1};
data_clean.time{1}	= b(selt);
b					= data_clean.trial{1};
b					= b(:,selt);
data_clean.trial{1} = b;
b					= data_clean.sync{1};
b					= b(selt);
data_clean.sync{1}	= b;

dt = [0 diff(selt)];
rmsample = find(dt>0);
data_clean.trigger{1}	= data_clean.trigger{1}-rmsample;

%%
if checkClean
	a				= [data_trigger.trigger{1}];
	t				= data_trigger.time{1}; % s
	selt			= t>a(1)/data_trigger.fsample-10 & t<a(end)/data_trigger.fsample+10; % within 10 s of on- or offset
	
	figure
	subplot(211)
	plot(t,data_raw.trial{1}(1,:));
	
	hold on
	plot(t(selt),data_raw.trial{1}(1,selt));
end

function data_down = resampledata(data_clean)

cfg				= [];
cfg.resamplefs	= 10;
cfg.demean		= 'yes';
data_down		= ft_resampledata(cfg, data_clean);
data_down.sync{1}		= transpose(resample(transpose(data_clean.sync{1}),cfg.resamplefs,data_clean.fsample));
data_down.trigger{1}	= round(data_clean.trigger{1}/data_clean.fsample*cfg.resamplefs);

function data_trans		= transformdata(data_nophys)

%% Transform
cfg = [];
cfg.target	= {'O2Hb', 'HHb'};
cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all?, all NIRS channels are selected then
data_trans	= ft_transform_ODs(cfg, data_nophys);

function data_epoched = epochdata(data_rcs,StmMtx,side)
%% Epoching
% t = getnirseventcode_AVCI;

cfg						= [];
cfg.trialdef.eventvalue = 1; % read conditions
cfg.trialdef.prestim    = 5; % in seconds
cfg.trialdef.poststim   = 20; % in seconds

trg	= data_rcs.trigger{1}; % trigger indices for onset of the standard
trg = trg-2*data_rcs.fsample;

t = (StmMtx(:,1)<0)+1;
sel = t==side;
trg = trg(sel);
trl						= NaN(length(trg)-2,4);

for stdidx = 1:length(trg)
	begsample		= trg(stdidx) - cfg.trialdef.prestim*data_rcs.fsample;
	endsample		= trg(stdidx) + cfg.trialdef.poststim*data_rcs.fsample - 1;
	offset			= -cfg.trialdef.prestim*data_rcs.fsample;
	trigger			= t(stdidx); % remember the trigger (=condition) for each trial
	trl(stdidx, :)	= [round([begsample endsample offset])  trigger];
end

cfg.trl			= trl;
cfg.lpfilter	= 'yes';                              % apply lowpass filter
cfg.lpfreq		= 0.08;                                 % lowpass
cfg.hpfilter	= 'yes';                              % apply highpass filter
cfg.hpfreq		= 0.008;                                 %
data_epoched	= ft_redefinetrial(cfg,data_rcs);


function data_rcs  = getdeep(data_trans)

data				= data_trans;
data	= nirs_optodedistance(data);

d					= data.opto.fiberdistance;

dat					= data.trial{:};
seldeep				= d>=3;
deep				= dat(seldeep,:);
time				= data.time{:};
data_rcs			= data;
data_rcs.time{1}	= time;
data_rcs.trial{1}	= deep;
label	= data.label; % transformed channel label

data_rcs.label		= label(seldeep);
