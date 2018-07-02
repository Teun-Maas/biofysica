function handles = trialShowMinor(handles)
% Modify to your own preference

% need for automatic calibration
% need for automatic detection


%% calibrated traces
% history in gray
% current black bold
handles.stim	= handles.trial(handles.cfg.trial).stim;
handles			= stimShowMinor(handles); % This should work ('Minor' not needed)
handles			= datShowMinor(handles); % probably needs adjustments

%% Linear Regression
% subplot(132)
% subplot(133)

%% Start drawing
drawnow;
