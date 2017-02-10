function handles = trialShow(handles)
% Modify to your own preference

% need for automatic calibration
% need for automatic detection


%% calibrated traces
% history in gray
% current black bold
handles.stim	= handles.trial(handles.cfg.trial).stim;
handles			= stimShow(handles);
handles			= datShow(handles);

%% Linear Regression
% subplot(132)
% subplot(133)

%% Start drawing
drawnow;
