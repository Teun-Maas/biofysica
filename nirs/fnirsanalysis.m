% basic NIRS analysis using fieldtrip functions
% requires the optodetemplates.xml in the added path
clear all


warning off

w			= what('LR-04-2015-06-17_active');
DataFolder	= w.path;
cd(DataFolder)
fname		= 'LR-04-2015-06-17-0001.oxy3';

%%
%read the raw data from file as one long continuous segment without any additional filtering
cfg				= [];
cfg.dataset      = fname;
data_raw = ft_preprocessing(cfg); 
save([DataFolder 'data_raw.mat'], 'data_raw');
load([DataFolder 'data_raw.mat']);

%%
% %plot selected channel and visualize event onset - raw data
figure(1)
chansel  = [5 6]; 
plot(data_raw.time{1}, data_raw.trial{1}(chansel(1), :),'Color',[1 0 1])
hold on
plot(data_raw.time{1}, data_raw.trial{1}(chansel(2), :))
xlabel('time (s)')
%ylabel('channel amplitude (uV)')
legend('[844nm]','[762nm]')
%legend(data_raw.label(chansel(1)),leg2)

%get indices for onset of triggered events
AD1 = data_raw.trial{1}(97,:); % standard
AD2 = data_raw.trial{1}(98,:); % deviant
n_sta		= [0 diff(AD1)];
sta		= find(n_sta>1.24); %standard 
nstim_sta	= numel(sta);
n_dev		= [0 diff(AD2)];
dev		= find(n_dev>1.24); %deviant 
nstim_dev	= numel(dev);
% draw as a line in the figure
pa_verline(data_raw.time{1}(sta));
drawnow
hold off

return
%%
%Segmenting continuous data into trials
% using my own function for trial selection (necessary due to triggering in AD channels)

%standard epoching

cfg = [];
cfg.dataset  = fname; 
cfg.trialdef.eventvalue = 1; % read conditions
cfg.trialdef.prestim    = 5; % in seconds
cfg.trialdef.poststim   = 15; % in seconds

trl = [];
for i=1:length(sta)-2
      begsample     = sta(i) - cfg.trialdef.prestim*data_raw.hdr.Fs;
      endsample     = sta(i) + cfg.trialdef.poststim*data_raw.hdr.Fs - 1;
      offset        = -cfg.trialdef.prestim*data_raw.hdr.Fs;  
      trigger       = 1; % remember the trigger (=condition) for each trial
      trl(end+1, :) = [round([begsample endsample offset])  trigger]; 
end

cfg.trl = trl;
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 10;                                 % lowpass
%cfg.hpfilter   = 'yes';                              % apply highpass filter
%cfg.hpfreq     = 0.02;                                 % 
data_sta_epoched = ft_preprocessing(cfg);

%deviant epoching
cfg = [];
%cfg.dataset  = 'LR-AR-2015-04-21-0001.oxy3'; 
cfg.dataset  = fname; 
cfg.trialdef.eventvalue = 2; % read conditions
cfg.trialdef.prestim    = 5; % in seconds
cfg.trialdef.poststim   = 15; % in seconds

trl = [];
for i=1:length(dev)
      begsample     = dev(i) - cfg.trialdef.prestim*data_raw.hdr.Fs;
      endsample     = dev(i) + cfg.trialdef.poststim*data_raw.hdr.Fs - 1;
      offset        = -cfg.trialdef.prestim*data_raw.hdr.Fs;  
      trigger       = 2; % remember the trigger (=condition) for each trial
      trl(end+1, :) = [round([begsample endsample offset])  trigger]; 
end
cfg.trl = trl;

cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 10;                                 % lowpass
%cfg.hpfilter   = 'yes';                              % apply highpass filter
%cfg.hpfreq     = 0.02;                                 % 
data_dev_epoched = ft_preprocessing(cfg);

%
figure(4)
chansel  = [1 2]; 
plot(data_sta_epoched.time{1}, data_sta_epoched.trial{1}(chansel(1), :),'Color',[1 0 1])
hold on
plot(data_sta_epoched.time{1}, data_sta_epoched.trial{1}(chansel(2), :))
xlabel('time (s)')
%ylabel('channel amplitude (uV)')
legend('[844nm]','[762nm]')
%legend(data_raw.label(chansel(1)),leg2)
hold off

figure(5)
chansel  = [1 2]; 
plot(data_dev_epoched.time{1}, data_dev_epoched.trial{1}(chansel(1), :),'Color',[1 0 1])
hold on
plot(data_dev_epoched.time{1}, data_dev_epoched.trial{1}(chansel(2), :))
xlabel('time (s)')
%ylabel('channel amplitude (uV)')
legend('[844nm]','[762nm]')
%legend(data_raw.label(chansel(1)),leg2)
hold off

save([DataFolder 'data_sta_epoched.mat'], 'data_sta_epoched','-v7.3');
save([DataFolder 'data_dev_epoched.mat'], 'data_dev_epoched');
load([DataFolder 'data_sta_epoched.mat']);
load([DataFolder 'data_dev_epoched.mat']);



%%
%convert data, i.e. changes in optical density to concentration changes
cfg = [];
cfg.target = {'O2Hb', 'HHb'};
cfg.channel = 'nirs'; % e.g. one channel incl. wildcards, you can also use ?all?, all NIRS channels are selected then
data_sta_trans = ft_transform_ODs(cfg, data_sta_epoched);
data_dev_trans = ft_transform_ODs(cfg, data_dev_epoched);

save([DataFolder 'data_sta_trans.mat'], 'data_sta_trans','-v7.3');
save([DataFolder 'data_dev_trans.mat'], 'data_dev_trans');

load([DataFolder 'data_sta_trans.mat']);
load([DataFolder 'data_dev_trans.mat']);


%%
%plot selected channel and visualize event onset
figure(6)
chansel  = [1 2]; 
plot(data_sta_trans.time{2}, data_sta_trans.trial{6}(chansel(1), :),'Color',[1 0 1])
hold on
plot(data_sta_trans.time{2}, data_sta_trans.trial{6}(chansel(2), :))
xlabel('time (s)')
%ylabel('channel amplitude (uV)')
legend('HbDataO2Hb','HbDataHHb')
hold off

figure(7)
chansel  = [1 2]; 
plot(data_dev_trans.time{2}, data_dev_trans.trial{6}(chansel(1), :),'Color',[1 0 1])
hold on
plot(data_dev_trans.time{2}, data_dev_trans.trial{6}(chansel(2), :))
xlabel('time (s)')
%ylabel('channel amplitude (uV)')
legend('HbDataO2Hb','HbDataHHb')
hold off

%%
% averaging
cfg = [];

%cfg.channel = 'gui'; % opens a gui to select channels
data_sta_avg = ft_timelockanalysis(cfg, data_sta_trans);
save([DataFolder 'data_sta_avg.mat'], 'data_sta_avg');
data_dev_avg = ft_timelockanalysis(cfg, data_dev_trans);
save([DataFolder 'data_dev_avg.mat'], 'data_dev_avg');


%to perform baseline normalization
cfg.baseline = [-5 0];
data_sta_avg_basel = ft_timelockbaseline(cfg, data_sta_avg);
save([DataFolder 'data_sta_avg_basel.mat'], 'data_sta_avg_basel');
data_dev_avg_basel = ft_timelockbaseline(cfg, data_dev_avg);
save([DataFolder 'data_dev_avg_basel.mat'], 'data_dev_avg_basel');
%save timelock
load([DataFolder 'data_sta_avg_basel.mat']);
load([DataFolder 'data_dev_avg_basel.mat']);


%%
% plotting
% %plot selected channel
figure(8)
chansel  = [45 46]; 
plot(data_sta_avg_basel.time,data_sta_avg_basel.avg(chansel(1), :),'r');
hold on
plot(data_sta_avg_basel.time,data_sta_avg_basel.avg(chansel(2), :),'b');
xlabel('time (s)')
ylabel('channel amplitude')
ylim([-0.2 0.2]);

figure(9)
chansel  = [45 46]; 
plot(data_dev_avg_basel.time,data_dev_avg_basel.avg(chansel(1), :),'r');
hold on
plot(data_dev_avg_basel.time,data_dev_avg_basel.avg(chansel(2), :),'b');
xlabel('time (s)')
ylabel('channel amplitude')
ylim([-0.2 0.2]);


%plot grids
figure(10) %left side
label = data_dev_avg_basel.label;
s={'Rx1a-Tx1 [O2Hb]', 'Rx2a-Tx1 [O2Hb]','Rx2a-Tx2 [O2Hb]','Rx1a-Tx3 [O2Hb]','Rx3a-Tx1 [O2Hb]','Rx2a-Tx4 [O2Hb]','Rx4a-Tx2 [O2Hb]','Rx3a-Tx3 [O2Hb]','Rx3a-Tx4 [O2Hb]','Rx4a-Tx4 [O2Hb]','Rx4b-Tx3 [O2Hb]','Rx3a-Tx5 [O2Hb]','Rx1b-Tx4 [O2Hb]','Rx4a-Tx6 [O2Hb]','Rx4b-Tx5 [O2Hb]','Rx1b-Tx5 [O2Hb]','Rx1b-Tx6 [O2Hb]','Rx4b-Tx7 [O2Hb]','Rx2b-Tx5 [O2Hb]','Rx1b-Tx8 [O2Hb]','Rx3b-Tx6 [O2Hb]','Rx2b-Tx7 [O2Hb]','Rx2b-Tx8 [O2Hb]','Rx3b-Tx8 [O2Hb]'};
s2={'Rx1a-Tx1 [HHb]', 'Rx2a-Tx1 [HHb]','Rx2a-Tx2 [HHb]','Rx1a-Tx3 [HHb]','Rx3a-Tx1 [HHb]','Rx2a-Tx4 [HHb]','Rx4a-Tx2 [HHb]','Rx3a-Tx3 [HHb]','Rx3a-Tx4 [HHb]','Rx4a-Tx4 [HHb]','Rx4b-Tx3 [HHb]','Rx3a-Tx5 [HHb]','Rx1b-Tx4 [HHb]','Rx4a-Tx6 [HHb]','Rx4b-Tx5 [HHb]','Rx1b-Tx5 [HHb]','Rx1b-Tx6 [HHb]','Rx4b-Tx7 [HHb]','Rx2b-Tx5 [HHb]','Rx1b-Tx8 [HHb]','Rx3b-Tx6 [HHb]','Rx2b-Tx7 [HHb]','Rx2b-Tx8 [HHb]','Rx3b-Tx8 [HHb]'};
s_all=[s; s2];
for i = 1:length(s)
    %plot selected channel
    h(i)=subplot(7,7,i*2);
    chan=find(strcmp(label,s(i))); % 
    plot(data_dev_avg_basel.time,data_dev_avg_basel.avg(chan, :),'r');
    hold on
    chan=find(strcmp(label,s2(i))); % 
    plot(data_dev_avg_basel.time,data_dev_avg_basel.avg(chan, :),'b');
    title(s{i}(1:8));
end
linkaxes(h)
ylim([-0.2 0.2]);
xlim([-5 15]);
xlabel('time (s)')
ylabel('channel amplitude')

figure(11) % right side
label = data_dev_avg_basel.label;
s={'Rx5a-Tx9 [O2Hb]', 'Rx6a-Tx9 [O2Hb]','Rx6a-Tx10 [O2Hb]','Rx5a-Tx11 [O2Hb]','Rx7a-Tx9 [O2Hb]','Rx6a-Tx12 [O2Hb]','Rx8a-Tx10 [O2Hb]','Rx7a-Tx11 [O2Hb]','Rx7a-Tx12 [O2Hb]','Rx8a-Tx12 [O2Hb]','Rx8b-Tx11 [O2Hb]','Rx7a-Tx13 [O2Hb]','Rx5b-Tx12 [O2Hb]','Rx8a-Tx14 [O2Hb]','Rx8b-Tx13 [O2Hb]','Rx5b-Tx13 [O2Hb]','Rx5b-Tx14 [O2Hb]','Rx8b-Tx15 [O2Hb]','Rx6b-Tx13 [O2Hb]','Rx5b-Tx16 [O2Hb]','Rx7b-Tx14 [O2Hb]','Rx6b-Tx15 [O2Hb]','Rx6b-Tx16 [O2Hb]','Rx7b-Tx16 [O2Hb]'};
s2={'Rx5a-Tx9 [HHb]', 'Rx6a-Tx9 [HHb]','Rx6a-Tx10 [HHb]','Rx5a-Tx11 [HHb]','Rx7a-Tx9 [HHb]','Rx6a-Tx12 [HHb]','Rx8a-Tx10 [HHb]','Rx7a-Tx11 [HHb]','Rx7a-Tx12 [HHb]','Rx8a-Tx12 [HHb]','Rx8b-Tx11 [HHb]','Rx7a-Tx13 [HHb]','Rx5b-Tx12 [HHb]','Rx8a-Tx14 [HHb]','Rx8b-Tx13 [HHb]','Rx5b-Tx13 [HHb]','Rx5b-Tx14 [HHb]','Rx8b-Tx15 [HHb]','Rx6b-Tx13 [HHb]','Rx5b-Tx16 [HHb]','Rx7b-Tx14 [HHb]','Rx6b-Tx15 [HHb]','Rx6b-Tx16 [HHb]','Rx7b-Tx16 [HHb]'};
s_all=[s; s2];
for i = 1:length(s)
    %plot selected channel
    h(i)=subplot(7,7,i*2);
    chan=find(strcmp(label,s(i))); % 
    plot(data_dev_avg_basel.time,data_dev_avg_basel.avg(chan, :),'r');
    hold on
    chan=find(strcmp(label,s2(i))); % 
    plot(data_dev_avg_basel.time,data_dev_avg_basel.avg(chan, :),'b');
    title(s{i}(1:8));
end
linkaxes(h)
ylim([-0.2 0.2]);
xlim([-5 15]);
xlabel('time (s)')
ylabel('channel amplitude')



%plot grid
figure(12) %left
label = data_sta_avg_basel.label;
s={'Rx1a-Tx1 [O2Hb]', 'Rx2a-Tx1 [O2Hb]','Rx2a-Tx2 [O2Hb]','Rx1a-Tx3 [O2Hb]','Rx3a-Tx1 [O2Hb]','Rx2a-Tx4 [O2Hb]','Rx4a-Tx2 [O2Hb]','Rx3a-Tx3 [O2Hb]','Rx3a-Tx4 [O2Hb]','Rx4a-Tx4 [O2Hb]','Rx4b-Tx3 [O2Hb]','Rx3a-Tx5 [O2Hb]','Rx1b-Tx4 [O2Hb]','Rx4a-Tx6 [O2Hb]','Rx4b-Tx5 [O2Hb]','Rx1b-Tx5 [O2Hb]','Rx1b-Tx6 [O2Hb]','Rx4b-Tx7 [O2Hb]','Rx2b-Tx5 [O2Hb]','Rx1b-Tx8 [O2Hb]','Rx3b-Tx6 [O2Hb]','Rx2b-Tx7 [O2Hb]','Rx2b-Tx8 [O2Hb]','Rx3b-Tx8 [O2Hb]'};
s2={'Rx1a-Tx1 [HHb]', 'Rx2a-Tx1 [HHb]','Rx2a-Tx2 [HHb]','Rx1a-Tx3 [HHb]','Rx3a-Tx1 [HHb]','Rx2a-Tx4 [HHb]','Rx4a-Tx2 [HHb]','Rx3a-Tx3 [HHb]','Rx3a-Tx4 [HHb]','Rx4a-Tx4 [HHb]','Rx4b-Tx3 [HHb]','Rx3a-Tx5 [HHb]','Rx1b-Tx4 [HHb]','Rx4a-Tx6 [HHb]','Rx4b-Tx5 [HHb]','Rx1b-Tx5 [HHb]','Rx1b-Tx6 [HHb]','Rx4b-Tx7 [HHb]','Rx2b-Tx5 [HHb]','Rx1b-Tx8 [HHb]','Rx3b-Tx6 [HHb]','Rx2b-Tx7 [HHb]','Rx2b-Tx8 [HHb]','Rx3b-Tx8 [HHb]'};
s_all=[s; s2];
for i = 1:length(s)
    %plot selected channel
    h(i)=subplot(7,7,i*2);
    chan=find(strcmp(label,s(i))); % 
    plot(data_sta_avg_basel.time,data_sta_avg_basel.avg(chan, :),'r');
    hold on
    chan=find(strcmp(label,s2(i))); % 
    plot(data_sta_avg_basel.time,data_sta_avg_basel.avg(chan, :),'b');
    title(s{i}(1:8));
end
linkaxes(h)
ylim([-0.2 0.2]);
xlim([-5 15]);
xlabel('time (s)')
ylabel('channel amplitude')


%plot grid
figure(13) % right
label = data_sta_avg_basel.label;
s={'Rx5a-Tx9 [O2Hb]', 'Rx6a-Tx9 [O2Hb]','Rx6a-Tx10 [O2Hb]','Rx5a-Tx11 [O2Hb]','Rx7a-Tx9 [O2Hb]','Rx6a-Tx12 [O2Hb]','Rx8a-Tx10 [O2Hb]','Rx7a-Tx11 [O2Hb]','Rx7a-Tx12 [O2Hb]','Rx8a-Tx12 [O2Hb]','Rx8b-Tx11 [O2Hb]','Rx7a-Tx13 [O2Hb]','Rx5b-Tx12 [O2Hb]','Rx8a-Tx14 [O2Hb]','Rx8b-Tx13 [O2Hb]','Rx5b-Tx13 [O2Hb]','Rx5b-Tx14 [O2Hb]','Rx8b-Tx15 [O2Hb]','Rx6b-Tx13 [O2Hb]','Rx5b-Tx16 [O2Hb]','Rx7b-Tx14 [O2Hb]','Rx6b-Tx15 [O2Hb]','Rx6b-Tx16 [O2Hb]','Rx7b-Tx16 [O2Hb]'};
s2={'Rx5a-Tx9 [HHb]', 'Rx6a-Tx9 [HHb]','Rx6a-Tx10 [HHb]','Rx5a-Tx11 [HHb]','Rx7a-Tx9 [HHb]','Rx6a-Tx12 [HHb]','Rx8a-Tx10 [HHb]','Rx7a-Tx11 [HHb]','Rx7a-Tx12 [HHb]','Rx8a-Tx12 [HHb]','Rx8b-Tx11 [HHb]','Rx7a-Tx13 [HHb]','Rx5b-Tx12 [HHb]','Rx8a-Tx14 [HHb]','Rx8b-Tx13 [HHb]','Rx5b-Tx13 [HHb]','Rx5b-Tx14 [HHb]','Rx8b-Tx15 [HHb]','Rx6b-Tx13 [HHb]','Rx5b-Tx16 [HHb]','Rx7b-Tx14 [HHb]','Rx6b-Tx15 [HHb]','Rx6b-Tx16 [HHb]','Rx7b-Tx16 [HHb]'};
s_all=[s; s2];
for i = 1:length(s)
    %plot selected channel
    h(i)=subplot(7,7,i*2);
    chan=find(strcmp(label,s(i))); % 
    plot(data_sta_avg_basel.time,data_sta_avg_basel.avg(chan, :),'r');
    hold on
    chan=find(strcmp(label,s2(i))); % 
    plot(data_sta_avg_basel.time,data_sta_avg_basel.avg(chan, :),'b');
    title(s{i}(1:8));
end
linkaxes(h)
ylim([-0.2 0.2]);
xlim([-5 15]);
xlabel('time (s)')
ylabel('channel amplitude')

saveas(figure(10), [DataFolder 'deviant_TP10Hz_left.fig']);
saveas(figure(11), [DataFolder 'deviant_TP10Hz_right.fig']);
saveas(figure(12), [DataFolder 'standard_TP10Hz_left.fig']);
saveas(figure(13), [DataFolder 'standard_TP10Hz_right.fig']);
