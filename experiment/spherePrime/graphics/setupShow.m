function handles = setupShow(handles)
% Modify to your own preference



% Trace: elevation vs azimuth
axes(handles.axes_xy); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
box off
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('azimuth (deg)');
ylabel('elevation (deg)');
title('2D trace');

% Stimulus-response plot+linear regression: azimuth
axes(handles.axes_regress_azimuth); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
unityline;
box off;
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('Target (deg)');
ylabel('Response (deg)');
title('Azimuth');

% Stimulus-response plot+linear regression: elevation
axes(handles.axes_regress_elevation); %#ok<*NASGU>
cla
hold on
axis([-120 120 -120 120]);
axis square;
unityline;
box off
set(gca,'XTick',-90:30:90,'YTick',-90:30:90,'TickDir','out');
xlabel('Target (deg)');
ylabel('Response (deg)');
title('Elevation');

axes(handles.axes_position); %#ok<*NASGU>
cla
hold on
box off
set(gca,'YTick',-90:30:90,'TickDir','out');
xlabel('Time (s)');
ylabel('Position (deg)');
ylim([-100 +100]);

axes(handles.axes_velocity); %#ok<*NASGU>
cla
hold on
box off
set(gca,'YTick',0:50:200,'TickDir','out');
xlabel('Time (s)');
ylabel('Velocity (deg)');
ylim([0 200]);
