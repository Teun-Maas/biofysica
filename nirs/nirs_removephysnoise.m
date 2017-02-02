function data = nirs_removephysnoise(data,varargin)
% DATA = NIRS_REMOVEPHYSNOISE(DATA)
%
% Remove physiological noise, such as heartbeat, respiration, and Mayer
% waves from NIRS DATA.
% Makes use of FieldTrip DATA structure and NAP's removeheartbeat function.
%
% See also: Nirs Analysis Package, Fieldtrip

%% Initialization
dispFlag = keyval('disp',varargin,false);


d = data.trial{1}';


d			= removeheartbeat(d,0.1); % heartbeat
d			= removeheartbeat(d,0.1,2); % respiration (0.2 Hz)
d			= removeheartbeat(d,0.1,[0.05 .2]); % Mayer Wave[0.05 .2]

data.trial{1} = d';
