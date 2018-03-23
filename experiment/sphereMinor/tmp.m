close all hidden
clearvars

%% Active X Control/Objects
% cfg.HF				= figure('Tag','ActXWin','name','ActiveX Window for TDT','numbertitle','off','menubar','none'); % Figure for ActiveX components
RZ6_1circuit = which('sphereMinor_RZ6_mp_velocity.rcx');

[zBus, err(1)]		= ZBUS(1); % zBus, number of racks
[RZ6_1, err(2)]	= RZ6(1,RZ6_1circuit); % Real-time acquisition


%% TDT status
RZ6_1Status	= RZ6_1.GetStatus;

%% Run zBus
zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).

%% Trigger event 1
% cfg.zBus.zBusTrigB(0, 0, 2); % start event 1, trial onset
zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;


pause(3)
%% Configuration

dur = 3; % s
nsamples = round(3*244.14)
% nsamples = 1000
RZ6_data		= NaN(nsamples,3);
tic;
for ii			= 1:3
    str = ['Data_' num2str(ii)]
    d = RZ6_1.ReadTagV(str,0,nsamples)';
    whos d
    RZ6_data(:,ii) = d;
    toc
end

%%

t = (1:nsamples)/244.14;
RZ6_data = RZ6_data*1700;

figure(2)
subplot(311)
plot(t,RZ6_data(:,1));
ylim([-10 10]);

subplot(312)
plot(t,RZ6_data(:,2));
ylim([-10 10]);

subplot(313)
plot(t,RZ6_data(:,3));
ylim([-10 10]);
