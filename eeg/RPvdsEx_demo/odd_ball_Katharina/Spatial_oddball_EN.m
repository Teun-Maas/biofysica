function Spatial_oddball
% Runs the spatial oddball experiment of Katherina Vogt

close all
clc
clearvars;

%-- Flags --%
Debug           =   1;  % Set to 1 to use the file "subj1/exp_subj1bl1.mat" for the experimental settings

%-- General Initialization --%
Path_and_name_rcx_file   = 'C:\Users\Elisabeth\Documents\MATLAB\odd_ball\Spatial_oddball_EN.rcx'; % The path + file name of the RPvdsEx circuit, for example the full path on the PC: 'C:\Users\Elisabeth\Documents\MATLAB\odd_ball\Spatial_oddball_EN.rcx'
Path_exp_settings = 'expblock'; % The folder where the file with the settings for each experiment is located. For example: C:\DATA\KV\spatial_plug\EXP\expblock\

if( Debug == 1 )
    %Folder		=	'DATA';
    Subj     	=	'1';
    Bl  	    =	'1';
elseif( Debug == 0 )
    %Folder		=	input('Enter the folder name: ','s');
    Subj    	=	input('Enter the subject code for the data file: ','s');
    Bl      	=	input('Enter the block number for the data file: ','s');
end

ExpName		=   [Path_exp_settings '\' 'subj' num2str(Subj) '\' 'exp_subj' num2str(Subj) 'bl' num2str(Bl) '.mat']; % The name of the file with the settings for this experiment: exp_subj1bl1.mat, the first number is the subject code, the second the block number.
load(ExpName);

%-- Intialization of variables --%
Par.RPpath		=	Path_and_name_rcx_file;
Par.Device		=	Sp.Grand_mat(:,7); % Device 1 or 3, this corresponds with
% PM2R multiplexer 0 cq 2 (1 or 3 is used as index
% of the array DeviceTable in line:
% RP.SetTagVal('Device',DeviceTable(Device));
Par.Bit			=	Sp.Grand_mat(:,6); % This determines together with the device number which speaker is activated
Par.Ndevice		=	length(Par.Device);
Par.Nbit		=	length(Par.Bit);
Ntrl            =   size(Sp.Grand_mat,1);
% Amp = 1; % amplitude of stimulus
ToneDur		=	Sp.ToneDur; %duration of stimulus
%Par.Npulse		=	Sp.Npulse; % number of pulses
Trig		=	Sp.Grand_mat(:,4); % The # of TTL pulses that needs to be generated. Normally this is for repeating
% a stimulus that is triggered by one external trigger (as the PulseTrain is used in the "Sound" tab of the circuit), 
% in the "Trigger" tab of oddball circuit it is used to save on the EEG
% recordings (the "digi" channel) which type of stimulus was presented 
% (the used speaker: 1 pulse corresponds to default: 0 degrees,
% 2 pulses indicate -30, etc, see Sp.Grand_mat column 4 and 5).
ISI 	    =	Sp.Grand_mat(:,8); % interstimulus time

%-- Setup TDT + loading circuit  --%
disp('Initializing...')
[RP,zBus,iniOK]	=	setuptdt_EN(Par.RPpath);
if iniOK % only continue when there are no errors
    disp('Initializing is ready!')
    
    %-- Main --%
    for k=1:Par.Ndevice  % Loop through the whole array with stimulus settings
        
        Device	=	Par.Device(k);  % Device and bit determine the speaker
        Bit	    =	Par.Bit(k);
        %ISI     =   Par.ISI(k);
        
        MUXSet_EN(RP,Device,Bit);	% Set the digital I/O of the RP6 to activate the correct bit on the correct device (PM2R multiplexer)
        RP.SoftTrg(1);			% Activate trigger 1 in the circuit: the stimulus with the settings shown below is generated
        disp(['Trial: ' num2str(k) ' of ' num2str(Ntrl),10, ' Device= ' num2str(Par.Device(k)) ' bit= ' num2str(Par.Bit(k)) ' ISI= ' num2str(ISI(k))])
        
        %RP.SetTagVal( 'Amp', Amp  ); % amplitude of stim
        RP.SetTagVal( 'HiTime', ToneDur); % duration of stim
        %RP.SetTagVal( 'NPulse', Npulse );
        RP.SetTagVal( 'Trig', Trig(k)); % Number of Pulses, input for the Pulse Train in the "Trigger" tab of the circuit
        RP.SetTagVal( 'LoTime', ISI(k)); % ISI
        
        zBus.zBusTrigA(0,0,10);	% zBusTrigA triggers several processor devices simultaneously.
        % The here used settings are: (trigger all devices, trig type = pulse, delay)
        % The delay before the trigger event occurs must be a minimum of 2msec
        % per rack.
        
        while( RP.GetTagVal('Active') )  % Wait while the stimulus is presented
            %-- hang during stimulus presentation --%
            pause(0.05);    % This is the time Matlab waits until checking
            % again if RP.GetTagVal('Active') is still true.
            % The tag "Active" is output of the PulseTrain in the
            % "Sound" tab of the circuit, 0 = waiting for trigger, 1 =
            % Output high, 2 = Output low, see the RPvdsEx manual.
            % The time specified in pause() only determines
            % the ISI when the time that RP.GetTagVal('Active')
            % needs to become false is shorter than this pause.
        end
        
        MUXSet_EN(RP,Device,0);  % Reset the three multiplexers. Remark: this will take a couple of msec. Problem?
    end
    disp('finished!')
end

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
        
        MUXClear_EN(RP); % Reset the three multiplexers. Question: is this necessary?
        pause(.05); % Some time for data transfer. Remark: this is only needed because of the SoftTrg in MUXClear?
    end
end

function MUXSet_EN(RP,Device,Channel)
% Set the digital I/O port  of the RP6 to activate the correct bit on the 
% correct device (PM2R multiplexer). If this function is called with the
% variable Channel filled, one of the speakers coupled to one of the three
% multiplexers is activated. If Channel is false all channels on the given
% Device are reset. (which devices will be addressed is determined by the
% settings for the digital I/O port in DeviceTable)
DeviceTable= [0 16 32];  % Original this was: [0 16 32 48], but there is no device at 48 in our TDT setup. The device at 16 is not in use.
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

    
function MUXClear_EN(RP)
% reset the three multiplexers (which devices will be addressed is determined
% by the settings for the digital I/O port of the RZ6 in DeviceTable)
DeviceTable= [0 16 32];
for Device=1:3
    RP.SetTagVal('Device',DeviceTable(Device));  % select the device
    RP.SetTagVal('SetReset',128);   % Give the tag in the circuit that is 
    % named SetReset the value needed for reset
    RP.SoftTrg(1); % Trigger soft trigger 1 in the circuit
    %	Note: see the above remark regarding SoftTrg
end

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
