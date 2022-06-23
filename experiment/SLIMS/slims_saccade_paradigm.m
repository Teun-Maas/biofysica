% SLIMS_SACCADE_PARADIGM
%
% Run the SACCADE experiment from the SLIMS project. For every experiment,
% the experimenter needs to change the parameters defining the filename (in
% the section Initialization).
%
% See also SLIMS_SACCADE_PREPROCESS

% To do:
% - GUI for filename
% - Separate alignment check (chin-rest)
% - LSL initialization in script
% - Saving data after each trial

%% Initialization
close all;
clearvars;

% GUI?
datadir			= 'C:\DATA\SS\JJ\Slims';
experimenter	= 'JJ';
subjectid		= '003';
sessionid		= '003';
formatout		= 'yy-mm-dd';
datetoday		= datestr(now,formatout);
fname			= fullfile([experimenter '-' subjectid '-' sessionid '-' datetoday '.mat']);
pcapture		= 'jesse_jansen'; 

%% Default psychtoolbox
% open a screen
Screen('Preference','SkipSyncTests',1);
PsychDefaultSetup(2); 
screens				= Screen('Screens');
screenNumber		= max(screens);
[window,windowRect] = PsychImaging('OpenWindow',screenNumber,0.2);
%search function PsychImaging
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
[xCenter, yCenter]	= RectCenter(windowRect); % center of screen
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
screenXdim			= 535; 
% dimensions of your screen
screenYdim			= 300;
distoscreen			= 0.5;
z					= distoscreen;
D					= z;
dotColor			= [1 0 0];
cdotColor			= [0 0 1];
dotSizePix			= 20;

%% Create target locations
step				= 0.2:0.4:1.8; % Locations are defined with respect to center [xCenter,yCenter] of screen, 0 = left,top edge, 2 = right,bottom edge
[stepx,stepy]		= meshgrid(step,step);
stepx				= stepx(:);
stepy				= stepy(:);
steps				= [stepx stepy];
ntrials				= length(steps);

%% Pseudo-randomize target presentation across trials
idx					= randperm(ntrials); % random indices
steps				= steps(idx,:); % implementing randomization

%% Randomize timing
fixtime				= randperm(500,ntrials); % 500 ms variation in additional fixation times, fixed time = 1s
tartime				= randperm(500,ntrials); % same for target  durations

% %% Start screen
% % Five dots are presented: 1 central and 4 edge targets. Subject is
% % required to fixate the central target. 
% baseRect			= [0 0 50 50];
% Rect				= CenterRectOnPointd(baseRect,1.9*xCenter,1.9*yCenter);
% 
% defscreen			= [0.2*xCenter 1.8*xCenter 0.2*xCenter 1.8*xCenter xCenter ; 1.8*yCenter 1.8*yCenter 0.2*yCenter 0.2*yCenter yCenter];
% Screen('DrawDots', window, defscreen, dotSizePix, dotColor, [], 2);
% Screen('Flip', window);
% defscreen			= [xCenter ; yCenter];
% WaitSecs(2);
% % rc = pupil_remote_control('131.174.140.178',50020);
% % r = rc.start_recording(pcapture);
% WaitSecs(0.1);
% 
% 
% pause; %  

%% LSL intialisation 
info_pl		= lsl_resolver('name=''Pupil Capture LSL Relay v2''');

l			= info_pl.list();
if isempty(l)
    error('no streams found');
end

for i		= 1:size(l,1)
    fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type);
end

if size(l,1)==1
    n		= 1;
else
    n		= input('enter pupil labs stream number to acquire: ');
end
plstr		= lsl_istream(info_pl{n});

metadata	= lsl_metadata_gaze(plstr);
%fprintf([metadata.as_xml() '\n']);
p			= metadata.as_struct(); 

info_ev		= lsl_resolver('type=''Digital Events @ lslder04'' and name=''Digital Events 1''');
l			= info_ev.list();
if isempty(l)
    error('no streams found');
end

for i = 1:size(l,1)
    fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type);
end

if size(l,1)==1
    n		= 1;
else
    n		= input('enter lslder stream number to acquire: ');
end
%      n=input('enter lslder stream number to acquire: ');
% n = 1;
evstr		= lsl_istream(info_ev{n});
ses			= lsl_session();
ses.add_stream(plstr);
ses.add_stream(evstr);
addlistener(plstr,'DataAvailable',@pl_listener);
addlistener(evstr,'DataAvailable',@ev_listener);

rc			= pupil_remote_control('131.174.140.178',50020);

ses.start();
WaitSecs(0.2);
rc.start_recording();
Screen('FillRect',window,[1,1,1],Rect);
Screen('Flip', window);

%% Run experiment!!!
% for every trial
t_fix		= NaN(ntrials,1);
t_start		= NaN(ntrials,1);
t_tar		= NaN(ntrials,1);
t_end		= NaN(ntrials,1);
targets		= NaN(ntrials,2); % target locations, referenced to center (see steps above)
for jj = 1:ntrials
	
	Screen('DrawDots', window,[xCenter 1.9*xCenter ; yCenter 1.9*yCenter], dotSizePix, [1,0,0], [], 2);
	Screen('FillRect',window,[0,0,0],Rect);
	t_fix(jj)	= Screen('Flip', window);
	t_start(jj) = GetSecs;
	WaitSecs(1+(fixtime(jj)/1000));
	
	Screen('FillRect',window,[1,1,1],Rect);
	Screen('Flip', window);
	WaitSecs(200/1000)
	
	Screen('DrawDots', window, [steps(jj,2)*xCenter  1.9*xCenter; steps(jj,1)*yCenter 1.9*yCenter], dotSizePix,[1,0,0], [], 2);
	Screen('FillRect',window,[0,0,0],Rect);
	t_tar(jj) = Screen('Flip', window);
	WaitSecs(1+(tartime(jj)/1000));
	t_end(jj) = GetSecs;
	
	Screen('FillRect',window,[1,1,1],Rect);
	Screen('Flip', window);
	WaitSecs(200/1000)
	
	targets(jj,:) = [steps(jj,2)*xCenter  steps(jj,1)*yCenter];
end
%%
rc.stop_recording();
pause(0.2);
ses.stop();
pldata	= plstr.read();
evdata	= evstr.read();
delete(ses);
delete(plstr);
delete(evstr);
delete(info_pl);
delete(info_ev);
delete(rc);
sca;

%% locations in angles 
% r       = rc.stop_recording;
% sel     = steps(:,1) == 1 & steps(:,2) ==1;
sel				= targets(:,1) == xCenter & targets(:,2) == yCenter;
xpixres			= screenXdim/screenXpixels;
ypixres			= screenYdim/screenYpixels;
xlen			= xpixres*targets(:,1);
ylen			= ypixres*targets(:,2);
xlen			= xlen - xlen(sel);
ylen			= ylen - ylen(sel);
xangle			= atand(xlen/(D*1000));
yangle			= atand(ylen/(D*1000));

%% Save
stimuli				= struct('target',[xangle yangle],'duration',[fixtime tartime],'norm_target',targets);
configuration		= struct('screendim',[screenXdim screenYdim],'screenpix',[screenXpixels screenYpixels],'distoscreen',D, 'targetinterval', [fixtime tartime],'pupildata',pldata,'eventdata',evdata,'metadata',metadata);
data				= struct('pupildata',pldata,'eventdata',evdata,'metadata',metadata);
filename			= fullfile(datadir,fname);
disp(filename);
save(filename,'stimuli','configuration','data');