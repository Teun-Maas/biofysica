 function AMtones_headphone_EN
% Produces left and right AM modulated sound for the ear phone.
% Author: Elisabeth Noordanus
%% Section: load settings
close all
clearvars;

%-- Flags --%
Testrun     = true;  % Set to true when doing a test run 

%-- General Initialization --%
Path_and_name_rcx_file   = 'C:\Users\Elisabeth\Documents\MATLAB\AMtones\2_AMtones_headphone_EN_runfrommatlab.rcx'; % The path + file name of the RPvdsEx circuit, for example the full path on the PC: 'C:\Users\Elisabeth\Documents\MATLAB\AMtones\2_AMtones_headphone_EN.rcx'
Path_exp_settings = 'expblock'; % Best to use the full pathname to prevent mixup with other dirs of the same name. The folder where the file with the settings for each experiment is located. For example: C:\DATA\KV\spatial_plug\EXP\expblock\
Saving_dir  = 'C:\Users\Elisabeth\Documents\MATLAB\';  % Address (see the line below) on the lab PC
% Saving_dir  = 'C:\Users\Elisabeth Noordanus\Documents\MATLAB\';  % End with backslash. Dir in which to save the produced file 
% Matlab always has write access to the MATLAB folder in C:\Users\ username \Documents\MATLAB\
% If you want to look at/change the file: double click it in Matlab,
% Windows has a problem with recognizing the correct file type if you
% didn't change the file association for the .mat extension.
if (~(7==exist(Saving_dir,'dir'))) % check if the folder exists
    % Saving the file with the executed experimental settings is only for
    % checking, that info should be identical to the file set below with
    % ExpName
    disp(['WARNING: The directory ' Saving_dir ' to save the experimental settings does not exist'])
end    
if (Testrun)
    Subj     	= 'test1';
    Bl  	    = 'test1';
else
    Subj    	= input('Enter the subject code for the data file: ','s');
    Bl      	= input('Enter the block number for the data file: ','s');
end

% Load the settings for eacht trial   
ExpName		=   [Path_exp_settings '\' 'AMtones-headphone_subj_' num2str(Subj) '_block_' num2str(Bl) '.mat']; % The name of the file with the settings for this experiment: exp_subj1bl1.mat, the first number is the subject code, the second the block number.
load(ExpName);
%-- Intialization of variables read from the file with settings per trial --%
Par.RPpath		= Path_and_name_rcx_file;
Par.Ntrials     = size(exp_settings.presented_stimuli,1);

%% Section: initialization  (but this section can't run stand alone because a function below is needed)
%-- Setup TDT + loading circuit  --%
disp('Initializing...')
[RP,zBus,iniOK]	=	setuptdt_EN(Par.RPpath);
if iniOK % check if the file with experimental settings exists
    if (~(2==exist(ExpName,'file'))) % check if the file exists
        disp(['ERROR: The folder ' ExpName ' with the experimental settings does not exist'])
        iniOK = false;
    end
end   
if iniOK % only continue when there are no errors
    disp('Initializing is ready!')
    %-- Set the variables that don't change during the experiment --%
    RP.SetTagVal( 'Att_A', exp_settings.Att_A); % Attenuation (db) of the left output
    RP.SetTagVal( 'Att_B', exp_settings.Att_B); % Attenuation (db) of the right output
    RP.SetTagVal( 'Stim_time_l', exp_settings.Stim_time_l); % Duration in msec when only the left sound is played
    RP.SetTagVal( 'Stim_time_r', exp_settings.Stim_time_r); % Duration in msec when only the right sound is played
    RP.SetTagVal( 'Stim_time_both', exp_settings.Stim_time_both); % Duration in msec when both sounds are played
    RP.SetTagVal( 'A_m_l', exp_settings.A_m_l); % Amplitude (V) of the modulation of the left tone
    RP.SetTagVal( 'A_m_r', exp_settings.A_m_r); % Amplitude (V) of the modulation of the right tone
    RP.SetTagVal( 'Ashift_m_l', exp_settings.Ashift_m_l); % Shift of the amplitude (V) of the modulation of the left tone
    RP.SetTagVal( 'Ashift_m_r', exp_settings.Ashift_m_r); % Shift of the amplitude (V) of the modulation of the right tone
    RP.SetTagVal( 'Ashift_c_l', exp_settings.Ashift_c_l); % Shift of the amplitude (V) of the carrier wave of the left tone
    RP.SetTagVal( 'Ashift_c_r', exp_settings.Ashift_c_r); % Shift of the amplitude (V) of the carrier wave of the right tone
    RP.SetTagVal( 'f_m_l', exp_settings.f_m_l); % Frequency (hz) of the modulation of the left tone
    RP.SetTagVal( 'f_m_r', exp_settings.f_m_r); % Frequency (hz) of the modulation of the right tone
    RP.SetTagVal( 'f_c_l', exp_settings.f_c_l); % Frequency (hz) of the carrier wave of the left tone
    RP.SetTagVal( 'f_c_r', exp_settings.f_c_r); % Frequency (hz) of the carrier wave of the right tone
    RP.SetTagVal( 'phase_m_l', exp_settings.phase_m_l); % Starting phase (rad) of the modulation of the left tone
    RP.SetTagVal( 'phase_m_r', exp_settings.phase_m_r); % Starting phase (rad) of the modulation of the right tone
    RP.SetTagVal( 'phase_c_l', exp_settings.phase_c_l); % Starting phase (rad) of the carrier wave of the left tone
    RP.SetTagVal( 'phase_c_r', exp_settings.phase_c_r); % Starting phase (rad) of the carrier wave of the right tone
    checkvalue = RP.GetTagVal('f_m_l');
    if (checkvalue)
        if (Testrun), display(['Modulation frequency left: ' num2str(checkvalue)]), end;
    else
        iniOK = false;
        disp('Error: values not correctly set');
    end
end
%% Section: Run experiment
presented_stimuli_check = zeros(Par.Ntrials,3);
if iniOK % only continue when there are no errors
    for k=1:Par.Ntrials  % Loop through all the trials
        ISI_1 = exp_settings.presented_stimuli(k,1); % Time (sec) to wait before presenting the sound which indicates on which side to concentrate
        ISI_2 = exp_settings.presented_stimuli(k,2); % Time (sec) to wait between presenting the sound which 
        concentrate_on_nr = exp_settings.presented_stimuli(k,3); % indicates the side to concentrate on and presenting the sound to both ears
        pause(ISI_1);                    
        % Pulse to indicate on which side the subject should concentrate
        concentrate_on = 'right';
        if (concentrate_on_nr == 1)
            RP.SoftTrg(1);			% Activate soft trigger 1 in the circuit: only the left stimulus will be heard for Stim_time_l seconds
            concentrate_on = 'left';
            ISI_2 = exp_settings.Stim_time_l/1000 + ISI_2;
        else
            RP.SoftTrg(2);			% Activate soft trigger 2 in the circuit: only the right stimulus will be heard for Stim_time_r seconds
            ISI_2 = exp_settings.Stim_time_r/1000 + ISI_2;
        end
        % MUXSet_EN(RP,Device,Bit);	% The setting of the multiplexers is done in the circuit. Set the digital I/O of the RP6 to activate the correct bit on the correct device (PM2R multiplexer)
       
        disp(['Trial: ' num2str(k) '. Concentrate on: ' concentrate_on '. ISI_1 was: ' num2str(ISI_1) '. ISI_2 is: ' num2str(ISI_2)])
        pause(ISI_2); 
        RP.SoftTrg(3);			% Activate soft trigger 3 in the circuit: both stimuli will be heard for Stim_time_both seconds
        % MUXSet_EN(RP,Device,0);  % Reset the three multiplexers. Not needed here. Remark: this will take a couple of msec. Problem?
        pause(exp_settings.Stim_time_both/1000);  % Wait while the sound is presented.
        presented_stimuli_check(k,:)=[ISI_1,ISI_2,concentrate_on_nr];
    end
    disp('finished!')
end
%% Section: Save settings
exp_settings_check = struct('Att_A',exp_settings.Att_A,'Att_B', exp_settings.Att_B,'Stim_time_l', exp_settings.Stim_time_l,'Stim_time_r', exp_settings.Stim_time_r, ...
    'Stim_time_both', exp_settings.Stim_time_both,'A_m_l', exp_settings.A_m_l,'A_m_r', exp_settings.A_m_r,'Ashift_m_l', exp_settings.Ashift_m_l,'Ashift_m_r', exp_settings.Ashift_m_r, ...
    'Ashift_c_l', exp_settings.Ashift_c_l,'Ashift_c_r', exp_settings.Ashift_c_r,'f_m_l', exp_settings.f_m_l,'f_m_r', exp_settings.f_m_r,'f_c_l', exp_settings.f_c_l, ...
    'f_c_r', exp_settings.f_c_r,'phase_m_l', exp_settings.phase_m_l,'phase_m_r', exp_settings.phase_m_r,'phase_c_l', exp_settings.phase_c_l,'phase_c_r', exp_settings.phase_c_r, ...
    'presented_stimuli', presented_stimuli_check);

% Save to file
FileName=[Saving_dir, 'AMtones-headphone-subj_',Subj,'_block_',Bl,'_',datestr(now, 'yyyymmdd_HHMMSS'),'.mat']; % add the date and time to the filename
save(FileName, 'exp_settings_check');

%-- Functions for interacting with the TDT in the EEG lab --%
function [RP,ZBus,iniOK] = setuptdt_EN(RPpath)
% This function creates activeX controls for the RZ6 and the zBus
iniOK = true;
%-- set up the ZBus --%
figure(100) % Figure for the activeX
set(gcf,'position',[1 1 1 1])
ZBus = actxcontrol('ZBUS.x',[1 1 .01 .01]); % The figure falls off the screen in the lower left corner with this setting. Try: [5 5 26 26] (see ActiveX manual), but what about the position setting in the line above?
if( ~ZBus.ConnectZBUS('GB') )
    iniOK = false;
    disp('Failed to init ZBus');
else
    %-- Set up RZ6 #1 via the Optical Gigabit interface to control MUXes (PM2Relays) --%
    %  The PM2Relays are 16 channel multiplexers that are used to switch the
    %  output signal from the RZ6 to one of maximal 16 speakers. Which of the
    %  multiplexers and speakers is activated depends on the setting of the
    %  digital I/O (25-pin serial cable that connects the RZ6 to the PM2Rs).
    RP	= actxcontrol('RPco.X',[1 1 .01 .01]);
    if( ~RP.ConnectRZ6('GB',1) ) % Connect to the RZ6, device number 1
        iniOK = false;
        disp('Failed to connect RP');
    else
        RP.Halt; % Stops any processing chains running on the RZ6 #1
        RP.ClearCOF; % Clears any circuit on the RZ6 #1
        disp(['Loading circuit ' RPpath]);
        RP.LoadCOF(RPpath); % Load the circuit file
        RP.Run; % Start the processor device's processing chain
        
        status=double(RP.GetStatus); % Gets the status
        if bitget(status,1)==0 % Checks for connection
            iniOK = false;
            disp('Error connecting to RZ6');
        elseif bitget(status,2)==0 % Checks for errors in loading circuit
            iniOK = false;
            disp('Error loading circuit');
        elseif bitget(status,3)==0 % Checks for errors in running circuit
            iniOK = false;
            disp('Error running circuit');
        else
            disp('Circuit loaded and running');
        end
        
        % MUXClear_EN(RP); % Reset the three multiplexers. Question: is this necessary?
        % pause(.05); % Some time for data transfer. Remark: this is only needed because of the SoftTrg in MUXClear?
    end
end


%{
% These functions are not used in this script
function MUXClear_EN(RP)
% reset the three multiplexers (which devices will be addressed is determined
% by the settings for the digital I/O port of the RZ6 in DeviceTable)
DeviceTable= [0 16 32];
for Device=1:3
    RP.SetTagVal('Device',DeviceTable(Device));  % select the device
    RP.SetTagVal('SetReset',128);   % Give the tag in the circuit that is 
    % named SetReset the value needed for reset
    RP.SoftTrg(1); % Trigger soft trigger 1 in the circuit. BE CAREFUL: change
    % this to the correct trigger! Must be implemented in the circuit. 
    %	Note: see the above remark regarding SoftTrg
end
function MUXSet_EN(RP,Device,Channel)
% Set the digital I/O port  of the RP6 to activate the correct bit on the 
% correct device (PM2R multiplexer). If this function is called with the
% variable Channel filled, one of the speakers coupled to one of the three
% multiplexers is activated. If Channel is false all channels on the given
% Device are reset. (which devices will be addressed is determined by the
% settings for the digital I/O port in DeviceTable)
DeviceTable= [0 16 32];  % Device 1, 2, 3
RP.SetTagVal('Device',DeviceTable(Device));  % set the device address
if Channel % activate a channel
    RP.SetTagVal('Chan',Channel-1); % select the channel, subtract 1 to obtain the correct bit since Grand_mat(:,6) counts from 1 and the bits are counted from 0
    RP.SetTagVal('SetReset',64);    % make it active
else % inactivate all channels on this device
    RP.SetTagVal('SetReset',128);
end
RP.SoftTrg(1); % Trigger soft trigger 1 in the circuit
    %	Note (see the ActiveX User Reference): Do not use software triggers for
    %	signal generation or acquisition that requires precise timing. Software
    %	triggers are affected by USB transfer times. Expect a 2-4 ms delay for
    %	each call to the processor device from the SoftTrg(). If multiple
    %	devices need to be triggered simultaneously use zBusTrigA/B().


function makedir_EN(DataPath)
% Create a new directory (not used currently in this script)
if( isunix )
    idx		=	strfind(DataPath,'/');
else
    idx		=	strfind(DataPath,'\');
end
NewPath	=	DataPath(1:idx(end)-1);
if( ~exist( NewPath,'dir' ) )
    mkdir(NewPath)
end
%}