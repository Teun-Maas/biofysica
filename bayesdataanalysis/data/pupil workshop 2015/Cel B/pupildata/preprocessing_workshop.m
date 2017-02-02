%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrofsen = 25;                      %number of presented sentences
exclsen = 5;                       %number of sentences excluded at start
exclpers = 15;                     %percentage blink for trial exclusion
timefreq = 50;                     %sample rate frequency
baseline_onset = 2*timefreq;       %start baseline at 2 sec
sentence_onset = 3*timefreq;       %start sentence at 3 sec
prompttone = 7*timefreq;           %start prompt tone at 7 sec 
                                   %shortest sentence (1.14 sec) + 3
taper = 3*timefreq;                %3 sec taper to avoid artifacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A time trace is created from -3 to 7 seconds to use for plotting
time = [-3:(1/timefreq):7-1/timefreq];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rawdata is transformed to a data set only containing the pupil traces
% (10 seconds long) of the trials of interest.   
[ldata,wdata] = size (rawdata);       %calculate dataset size
nroftr = wdata/5;                     %calculate nr of recorded traces
excltr = nroftr-(nrofsen-exclsen);    %exclude first traces
nroftr = nroftr-excltr;               %calculate nr of trails
ptrace = rawdata(1:prompttone+taper,(excltr*5)+3:5:wdata); % left pupil 
validtrial(1:nroftr) = 1;    % Valid trial list (1=included 0=excluded)                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the unprocessed traces    
figure
hold on
plot (time, ptrace)
xlim ([-3.0 7.0]);
hold off
MT = mean(ptrace');
figure
hold on
plot (time, MT)
xlim ([-3.0 7.0]);
hold off                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                     
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear interpolation of blinks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        %if a trace ends with a zero no interpolation is performed
        if isempty(ptrace(eblink,k)), break, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %linear interpolation 
            y = [ptrace(bblink,k) ptrace(eblink,k)]; % y difference
            gap = eblink - bblink;                   % x difference
            iss = (y(2)-y(1))/gap; %calculate the steps of interpolation
            yi = y(1);
            for m = 1: gap
                yi = [yi y(1)+m*iss];   % interpolation 
            end
            q=1;
            for h = bblink:eblink
                ptrace(h,k) = yi(q);    % implement interpolation
                q=q+1;
            end                        
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the interpolated traces 
figure
hold on
plot (time, ptrace)
xlim ([-3.0 7.0]);
hold off
MT = mean(ptrace');
figure
hold on
plot (time, MT)
xlim ([-3.0 7.0]);
hold off                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5-point moving average smoothing filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i =  1:nroftr
    trace = [];
    for j = 1:prompttone+taper
    trace = [trace ptrace(j,i)]; %#ok<AGROW>
    end
    for k = 3:prompttone-2+taper
      ptrace(k,i) = mean(trace(k-2:k+2));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the smoothed traces                           
figure
hold on
plot (time, ptrace)
xlim ([-3.0 7.0]);
hold off
MT = mean(ptrace');
figure
hold on
plot (time, MT)
xlim ([-3.0 7.0]);
hold off                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline correction and averaging to sentence onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%baseline correction
allbaseline = [];
for k= 1:nroftr;
    baseline = mean(ptrace(baseline_onset:sentence_onset,k));
        if validtrial(k) == 1
            allbaseline = [allbaseline baseline]; %#ok<AGROW>
        end
    for j = 1:prompttone+taper
    ptrace(j,k) = ptrace(j,k)- baseline; 
    end
end
meanbaseline = mean(allbaseline);
 
%avarage data of valid trials
avg_data = 0;
for i = 1:nroftr;
    if validtrial(i) == 1
        avg_data = ptrace(:,i)+ avg_data;
    end
end
avg_data = avg_data/length(find(validtrial==1));            %average trace
MPD = mean(avg_data(sentence_onset:prompttone));            %mean dilation
[PPD, latency] = max(avg_data(sentence_onset:prompttone));  %peak dilation
latency = latency/timefreq;                                 %peak latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the baseline corrected traces
figure
hold on
plot (time, ptrace)
xlim ([-3.0 7.0]);
hold off
figure
hold on
plot (time, avg_data)
xlim ([-1.0 5.0]);
ylim ([-0.05 0.5]);
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_data = {};
output_data{1,1} = time';               %time scale
output_data{1,2} = avg_data';           %mean trace
output_data{1,3} = ptrace;              %baseline corrected traces
output_data{1,4} = validtrial;          %list trials included
output_data{1,5} = percblinktrials;     %percentage of trials excluded
output_data{1,6} = allbaseline;         %all individual baselines
output_data{1,7} = meanbaseline;        %mean baseline
output_data{1,8} = MPD;                 %mean pupil dilation
output_data{1,9} = PPD;                 %peak pupil dilation
output_data{1,10} = latency;            %peak pupil latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%