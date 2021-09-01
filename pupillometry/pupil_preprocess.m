
function [time,data,Fs] = pupil_preprocess(time,data,Ti,varargin)
% In this script blinks are removed and traces are smoothed. Traces
% containing to many blinks are excluded. The output provides an average
% trace that is baseline corrected.
% The following processing steps are performed.
% PREPROCESSING  Collects pupil traces from the raw data set and
%                creates a list that codes for valid and invalid traces
%                based on the percentage of blinks.
% LINEAR INT     Fills in gaps created by blinks.
%                Note that blinks will also create gaps in eye-movement
%                data. This is not covered in this workshop.
% MOVING AVE     Moving average filter for smoothing data.
% BASELINE COR   Baseline correction is preferment for each trace.
% AVERAGING      Average over all traces is calculated.


 
rawFlag         = keyval('showRaw',varargin,false);
interpFlag      = keyval('showInterp',varargin,false);
smoothFlag      = keyval('showSmooth',varargin,false);
baselineFlag    = keyval('showBaseline',varargin,false);
avgFlag         = keyval('showAvg',varargin,false);
timeFlag         = keyval('showTime',varargin,false);
 
dispFlag         = keyval('disp',varargin,false);
if ~dispFlag
    rawFlag         = false;
    interpFlag      = false;
    smoothFlag      = false;
     
    timeFlag        = false;
end
 
%% Data
 
% [~,ntrials]       = size(data);           % Calculate dataset size
 
 
%% Timing
t   = time; % sample time (s)
dt  = diff(t); % time between samples (s)
f   = 1./dt; % sample rate (Hz)
 
% Fs = cellfun(@(x) x(1:nsamples),data,'UniformOutput',false)';
if timeFlag
    figure;
    subplot(211)
    plot(dt*1000,'ko-','MarkerFaceColor','w')
    ylabel('\delta t (ms)');
    xlabel('sample');
     
    subplot(223)
    x = 0.8/120:1/12000:(1.5/120);
    hist(dt,x)
    xlabel('\delta t (s)');
    ylabel('N');
     
    subplot(224)
    x = 0:1:240;
    hist(f,x)
    xlabel('sample rate (Hz)');
    ylabel('N');
end
%% Hard coded
Fs                      = mean(f); % ~119 Hz
 
%% Resample the data
Fs = 120;

[data,time] = resample(data,time,Fs);

%% Preprocessing
 
 
% data                  = data/10; % from pixels to mm, see pupillab manual/webpage/github wiki?
data = data';
[ntrials,nsamples]      = size(data) ;          % Calculate dataset size
 
 
%% Selection and time
 
% The data is transformed to a data set only containing the pupil traces
% (~9 seconds long) of the trials of interest.
 
if rawFlag, plotdata(data,time,101), title('Raw data');
end
 
%% Blink interpolation
% Linear interpolation of blinks
data = blinkinterp(data,time,ntrials);
 
 
%% Spike removal
% if exist('movart2clean','file')
%   data = movart2clean(data);
% end
%
if interpFlag, plotdata(data,time,102), title('Blinks removed');
end
% because not all blinks are removed
 
%% Bandpass & smooth
% 5-point moving average smoothing filter
% remove noise
% data        = movavg(data,5);
 
for trlIdx = 1:ntrials
    data(trlIdx,:) = highpass(data(trlIdx,:),'Fc',0.1,'Fs',Fs)'; % are these good settings??????? 0.1 Hz just means baseline removal
    data(trlIdx,:) = lowpass(data(trlIdx,:),'Fc',20,'Fs',Fs)'; % this might be too strict bc response itself is slow. bc takes 4 sec for pupil to dilate (impulse response). 
% 0.25 Hz (bc 4 s for max) so the data is slow. might set highpass to 0.01
% Hz. something further away. you don't want to remove the data with this
% filtering.
end
 
 
if smoothFlag 
    plotdata(data,time,103), title('Smoothed');
end

data = interp1(time,data,Ti,'linear',NaN);
time = Ti;



function data = blinkinterp(data,time,ntrials)
 
 
%%
for trlIdx = 1:ntrials
    trace               = data(trlIdx,:);
    sel                 = double(trace==0);
    sel                 = double(movavg(sel,8)>0);
    trace               = interp1(time(~sel),trace(~sel),time,'linear','extrap');
    data(trlIdx,:)      = trace;
end


 
function plotdata(data,time,fig)
figure(fig)
hold on
plot (time,data)
axis square;
box off
xlabel('Time (s)');
ylabel('Pupil diameter (px)');


function [data_out,data_in] = pupil_preprocess2(data_in,varargin)
% In this script blinks are removed and traces are smoothed. Traces
% containing to many blinks are excluded. The output provides an average
% trace that is baseline corrected.
% The following processing steps are performed.
% PREPROCESSING  Collects pupil traces from the raw data set and
%                creates a list that codes for valid and invalid traces
%                based on the percentage of blinks.
% LINEAR INT     Fills in gaps created by blinks.
%                Note that blinks will also create gaps in eye-movement
%                data. This is not covered in this workshop.
% MOVING AVE     Moving average filter for smoothing data.
% BASELINE COR   Baseline correction is preferment for each trace.
% AVERAGING      Average over all traces is calculated.


%% Initialization
if nargin<1
	% 	fname = '/Users/marcw/DATA/Sebastian Ausili/COMP_PUPILA_SEBA_tr10_rb_s0n0_IFFM__18Oct2016_1655.mat';
	
	fname = '/Users/marcw/DATA/Sebastian Ausili/COMP_Sancho_Pupil_tr30_rb_s0n0_IFFM__27Oct2016_1642.mat';
	load(fname);
	data_in = rec;
	close all;
end
if ischar(data_in)
	load(data_in);
	data_in = rec;
end
% Graphics Flags

rawFlag         = keyval('showRaw',varargin,true);
interpFlag		= keyval('showInterp',varargin,true);
smoothFlag      = keyval('showSmooth',varargin,true);
baselineFlag    = keyval('showBaseline',varargin,true);
avgFlag         = keyval('showAvg',varargin,true);
timeFlag         = keyval('showTime',varargin,true);

dispFlag         = keyval('disp',varargin,false);
if ~dispFlag
	rawFlag			= false;
	interpFlag		= false;
	smoothFlag		= false;
	baselineFlag	= false;
	avgFlag			= false;
	timeFlag		= false;
end

exclsen			= keyval('sentenceXclude',varargin,0);   % Sentence exclusion
exclperc		= keyval('blinkXclude',varargin,15);  % Percentage blink for trial exclusion

%% Data
data			= data_in.pupildata.diameter;
time			= data_in.pupildata.timestamp;
stim			= data_in.env_speech; % envelope

% [~,ntrials]		= size(data);           % Calculate dataset size


%% Timing
t	= time{1}; % sample time (s)
dt	= diff(t); % time between samples (s)
f	= 1./dt; % sample rate (Hz)
Fs	= mean(f); % mean sample frequency

% Fs = cellfun(@(x) x(1:nsamples),data,'UniformOutput',false)';
if timeFlag
	subplot(211)
	plot(dt*1000,'ko-','MarkerFaceColor','w')
	ylabel('\delta t (ms)');
	xlabel('sample');
	
	subplot(223)
	x = 0.8/120:1/12000:(1.5/120);
	hist(dt,x)
	xlabel('\delta t (s)');
	ylabel('N');
	
	subplot(224)
	x = 0:1:240;
	hist(f,x)
	xlabel('sample rate (Hz)');
	ylabel('N');
end
%% Hard coded
Fs						= 120;
Fsaudio					= 44100; % ?
baseline_onset			= 2*Fs; % Start baseline at 2 sec
sentence_onset_s		= 3; % Start sentence at 3 sec

sentence_onset			= 3*Fs; % Start sentence at 3 sec
taper					= round(2.5*Fs); % 3 sec taper to avoid artifacts
prompttone				= 6*Fs; % Start prompt tone at 6 sec `?????????


%% Preprocessing
% for ii = 1:ntrials
% nsamples(ii) = numel(data{ii});
% end
nsamples				= cellfun(@numel,data);
nsamples				= min(nsamples);
data					= cellfun(@(x) x(1:nsamples),data,'UniformOutput',false)';
data					= cell2mat(data)';

data					= data/10; % from ?? to mm, see pupillab manual/webpage/github wiki?
[nsamples,ntrials]		= size(data);           % Calculate dataset size


%% Selection and time
time		= 1:nsamples;
time		= time/Fs-sentence_onset_s;

% The data is transformed to a data set only containing the pupil traces
% (~9 seconds long) of the trials of interest.
excltr					= ntrials-(ntrials-exclsen);    %exclude first traces
ntrials					= ntrials-excltr;               %calculate nr of trails
idx						= 1:prompttone+taper;
data					= data(idx,excltr+1:ntrials); % left pupil
time					= time(idx);
validtrial(1:ntrials)	= 1;    % Valid trial list (1=included 0=excluded)

if rawFlag, plotdata(data,time,1), title('Raw data');
end

% [nsamples,ntrials]		= size(data);           % Calculate dataset size

%% Blink removal
% Calculate the percentage of blink. Blinks are shown as zero values.
for j = 1:ntrials
	blinkcount = 0;
	for i = baseline_onset:prompttone % from start baseline to prompt tone
		if data(i,j) == 0
			blinkcount = blinkcount+1;
		end
	end
	%if the blink count is higher than exclpers, then make the trial invalid
	if blinkcount/(prompttone-baseline_onset) > exclperc/100
		validtrial(j) = 0;
	end
end
percblinktrials = (1-(sum(validtrial)/ntrials))*100;

%% Blink interpolation
% Linear interpolation of blinks
data = blinkinterp(data,time,ntrials);


%% Spike removal
if exist('movart2clean','file')
	data = movart2clean(data);
end

if interpFlag, plotdata(data,time,2), title('Blinks removed');
end

%% Bandpass & smooth
% 5-point moving average smoothing filter
data		= movavg(data,5);

for trlIdx = 1:ntrials
	data(:,trlIdx) = highpass(data(:,trlIdx),'Fc',0.1,'Fs',Fs)';
	data(:,trlIdx) = lowpass(data(:,trlIdx),'Fc',2,'Fs',Fs)';
end


if smoothFlag, plotdata(data,time,3), title('Smoothed');
end

% %% Detrend
% for trlIdx = 1:ntrials
% 	y			= data(:,trlIdx);
% 	p			= polyfit(time,y',2);
% 	a			= y'-polyval(p,time);
% 	data(:,trlIdx) = a';
% end

%% Baseline correction and averaging to sentence onset
baseline		= mean(data(baseline_onset:sentence_onset,:));
data			= bsxfun(@minus,data,baseline);
if baselineFlag, plotdata(data,time,4), title('Baseline corrected');
end


%% Average
validtrial	= logical(validtrial);
avg_data	= nanmean(data(:,validtrial),2);
sd_data		= nanstd(data(:,validtrial),[],2);
se_data		= sd_data./sqrt(sum(validtrial)); % standard error



if avgFlag
	figure(6)
	hold on
	errorpatch(time,avg_data,se_data,[0.7 0 0]);
	% 	ylim ([-0.1 0.4]);
	xlim ([-1.5 6.0]);
	horline(0,'k-');
	verline([-1 0 5],'k:');
	axis square;
	xlabel('Time (s)');
	ylabel('Pupil diameter (mm)');
	box off;
	set(gca,...
		'YTick',0:0.1:0.3,...
		'XTick',-1:1:5,...
		'TickDir','out');
	
	
	%% HDR
	X		= zeros(size(time));
	sel		= time>=0& time <5;
	X(sel)	= 1;
	x		= nirs_hdrfunction(X,'Fs',Fs,'disp',false);
	y		= avg_data;
	b		= regstats(y,x,'linear','beta');
	y		= b.beta(2)*nirs_hdrfunction(X,'Fs',Fs,'disp',false);
	plot(time,y,'k-','LineWidth',3);
	plot(time,X);
	% 	savegraph(mfilename,'png');
	
	
end

%% Singular value decomposition
% close all
[u,s,v] = svd(detrend(data','linear'));
ev1		= u(1)*s(1)*v(:,1);
if avgFlag
	plot(time,ev1)
	hold on
end


%% Parameters & Output data
meanbaseline	= mean(baseline);
MPD				= mean(avg_data(sentence_onset:length(data)));% mean dilation
[PPD, latency]	= max(avg_data(sentence_onset:length(data))); % Peak pupil dilation
latency			= latency/Fs;                                 %peak latency

data_out.time				= time;                %time scale
data_out.avg				= avg_data';           %mean trace
data_out.se					= se_data';           %mean trace
data_out.data				= data;              %baseline corrected traces
data_out.valid				= validtrial;          %list trials included
data_out.percblinktrials	= percblinktrials;     %percentage of trials excluded
data_out.baseline			= baseline;         %all individual baselines
data_out.meanbaseline		= meanbaseline;        %mean baseline
data_out.MPD				= MPD;                 %mean pupil dilation
data_out.PPD				= PPD;                 %peak pupil dilation
data_out.latency			= latency;            %peak pupil latency
data_out.ev1				= ev1;


function data = blinkinterp(data,time,ntrials)
for trlIdx = 1:ntrials
	trace				= data(:,trlIdx);
	sel					= double(trace==0);
	sel					= double(movavg(sel,8)>0);
	trace				= interp1(time(~sel),trace(~sel),time,'linear','extrap');
	data(:,trlIdx)		= trace;
end

function plotdata(data,time,fig)
figure(fig)
hold on
plot (time,data)
axis square;
box off
xlabel('Time (s)');
ylabel('Pupil diameter (pixels?)');