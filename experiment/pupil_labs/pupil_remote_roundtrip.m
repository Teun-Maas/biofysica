% This example demonstrates how to send simple messages to the Pupil Remote plugin
%         'R' start recording with auto generated session name
%         'R rec_name' start recording and name new session name: rec_name
%         'r' stop recording
%         'C' start currently selected calibration
%         'c' stop currently selected calibration
%         'T 1234.56' Timesync: make timestamps count form 1234.56 from now on.
%         't' get pupil timestamp

% tested with jeromq and java 1.8
% GW/20170613-1

% do this before calling this script:
% MacOS
%javaclasspath('/Users/gunter/Documents/MATLAB/java/jeromq.jar');
% Ubuntu Linux
javaclasspath('/usr/share/java/jeromq.jar');

%rc = pupil_remote_control('dcn-eyebrain.local');
rc = pupil_remote_control('pupil-desktop.local');

r=rc.time_sync(0.0);

ntimes=1000;
rt=zeros(1,ntimes);
for i=1:ntimes
    tstart=tic;
    timestamp=rc.get_time_stamp();
    roundtrip=toc(tstart)*1000;
    rt(i)=roundtrip;
  %  fprintf(1,"pupil time stamp = %18.15f roundtrip time=%4.1f ms\n",timestamp, roundtrip);
end

r=rc.stop_recording;

delete(rc);


figure(gcf);
% probability density function
%histogram(rt,'Normalization','probability','BinWidth',0.1,'BinLimits',[0,10]);
% cumulative density function
histogram(rt,'Normalization','cdf','BinWidth',0.1,'BinLimits',[0,10]);

min_rt=min(rt)
max_rt=max(rt)
grid('on');
