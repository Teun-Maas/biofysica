function exLslPlotData2(evdata,pldata)
    
cts_evdata=lsl_correct_lsl_timestamps(evdata);
cts_pldata=lsl_correct_pupil_timestamps(pldata);
t0=min(cts_pldata(1),cts_evdata(1));
ev_y=evdata.Data;

%     Pupil Primitive data streams consist of 5 channels (lsl.cf_double64):
%         - diameter (-1.0 for gaze streams)
%         - confidence
%         - timestamp
%         - norm_pos.x
%         - norm_pos.y

pl_confidence=pldata.Data(2,:); % confidence

figure(101);
clf;
plot(cts_evdata-t0,ev_y,'o');
hold on;
plot(cts_pldata-t0,pl_confidence,'.');
set(gca,'XLim',[0 5]);
grid on;
grid minor;

end
