function data = pupilpreprocessing(rawdata,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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
% This MATLAB script is written by Thomas Koelewijn & Yang Wang 2015 for:
% The First Summer Workshop on The Application of Pupillometry in Hearing
% Science to assess Listening Effort VU University medical center,
% Amsterdam, The Netherlands, 5-6 August 2015

%% Initialization
if nargin<1
	fname = 'p14_celC.mat';
	load(fname)
end
freq	= keyval('freq',varargin,60);
% close all

nrofsen			= 25;                      %number of presented sentences
exclsen			= 5;                       %number of sentences excluded at start
exclpers		= 15;                     %percentage blink for trial exclusion
timefreq		= freq;                     %sample rate frequency
baseline_onset	= 2*timefreq;       %start baseline at 2 sec
sentence_onset	= 3*timefreq;       %start sentence at 3 sec
prompttone		= 7*timefreq;           %start prompt tone at 7 sec
%shortest sentence (1.14 sec) + 3
taper			= 3*timefreq;                %3 sec taper to avoid artifacts

% A time trace is created from -3 to 7 seconds to use for plotting
time			= -3:(1/timefreq):7-1/timefreq;

f = gcf;

%% Preprocessing
% The rawdata is transformed to a data set only containing the pupil traces
% (10 seconds long) of the trials of interest.
[~,wdata]	= size (rawdata);       %calculate dataset size
nroftr		= wdata/5;                     %calculate nr of recorded traces
excltr		= nroftr-(nrofsen-exclsen);    %exclude first traces
nroftr		= nroftr-excltr;               %calculate nr of trails
ptrace		= rawdata(1:prompttone+taper,(excltr*5)+3:5:wdata); % left pupil
validtrial(1:nroftr) = 1;    % Valid trial list (1=included 0=excluded)

% Plot the unprocessed traces
figure(f); subplot(221);
plotgraph(time,ptrace);
title('unprocessed traces');

% output_data= ptrace;
% 
% return
%% Linear interpolation of blinks

% Calculate the percentage of blink. Blinks are shown as zero values.
for j = 1:nroftr;
	blinkcount = 0;
	for i = baseline_onset:prompttone; % from start baseline to prompt tone
		if ptrace(i,j) == 0
			blinkcount = blinkcount+1;
		end
	end
	%if the blink count is higher than exclpers, then make the trial invalid
	if blinkcount/(prompttone-baseline_onset) > exclpers/100
		validtrial(j) = 0;
	end
end
percblinktrials = (1-(sum(validtrial)/nroftr))*100;
validtrial		= logical(validtrial);
data.validtrial = validtrial;

% Perform linear interpolation. This is done on all the traces.
% Traces with to many blinks can be excluded later.
for k = 1:nroftr; %For all traces
	for i= 6:prompttone+taper %For beginning to end of trace
		if ptrace(i,k) == 0 %If value is zero (blink)
			x = 0; %The counter x is set to zero
			bblink = i-5;  %Begin blink is set 5 samples before the blink
			%Defining the end of the blink
			for j = i+1:prompttone+taper %Start from next sample
				if ptrace(j,k) > 0          %If next sample is no blink
					x = x+1;                %x+1
				elseif ptrace(j,k) == 0     %if next sample is a blink
					x = 0;                  %x counter is reset
				end
				%x reaches 8 when for 8 samples in a row no zero value is found
				if x == 8
					eblink = j; %eblink is 8 samples after last zero value
					break
				end
			end
			
			%if a trace ends with a zero no interpolation is performed
			if isempty(ptrace(eblink,k)), break, end
			
			%linear interpolation
			y = [ptrace(bblink,k) ptrace(eblink,k)]; % y difference
			gap = eblink - bblink;                   % x difference
			iss = (y(2)-y(1))/gap; %calculate the steps of interpolation
			yi = y(1);
			for m = 1: gap
				yi = [yi y(1)+m*iss];   %#ok<AGROW> % interpolation
			end
			q=1;
			for h = bblink:eblink
				ptrace(h,k) = yi(q);    % implement interpolation
				q=q+1;
			end
		end
	end
end

% Plot the interpolated traces
figure(f); subplot(222);
plotgraph(time,ptrace);
title('interpolated traces');


%% 5-point moving average smoothing filter
ptrace = movavg(ptrace,5);

%Plot the smoothed traces
figure(f); subplot(223);
plotgraph(time,ptrace);
title('smoothed traces');



data.smoothed = ptrace;



%% baseline correction

baseline		= mean(ptrace(baseline_onset:sentence_onset,:));
ptrace			= bsxfun(@minus,ptrace,baseline);
meanbaseline	= mean(baseline(:,validtrial));

data.baseline = baseline;
data.baselinecorrected = ptrace;

% Plot the smoothed traces
figure(f); subplot(224);
plotgraph(time,ptrace);
title('baseline removal');
ylim([-1 1]);

%% average data of valid trials
avg_data				= mean(ptrace(:,validtrial),2);
MPD						= mean(avg_data(sentence_onset:prompttone));            %mean dilation
[PPD, latency]			= max(avg_data(sentence_onset:prompttone));  % peak dilation
latency					= latency/timefreq;                                 %peak latency

%% Output data
data.time				= time;                %time scale
data.average			= avg_data';           %mean trace
data.percblinktrials	= percblinktrials;     %percentage of trials excluded
data.MPD				= MPD;                 %mean pupil dilation
data.PPD				= PPD;                 %peak pupil dilation
data.latency			= latency;            %peak pupil latency
data.skewness			= skewness(data.smoothed(:));

function plotgraph(time,ptrace)
MT = mean(ptrace,2);

hold on
% for ii = 1:size(ptrace,2)
% 	plot (time, ptrace(:,ii)+ii*0.2,'k-','Color',[.7 .7 .7])
% 	hold on
% end
	plot (time, ptrace,'k-','Color',[.7 .7 .7])

plot (time, MT,'k-')

set(gca,'TickDir','out');
xlabel('Time (s)');
ylabel('Pupil size (mm)');
box off
xlim ([-3.0 7.0]);
title('unprocessed traces');

