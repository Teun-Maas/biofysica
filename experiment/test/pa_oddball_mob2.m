function Q = pa_oddball_mob2
% NIRS passive oddball experiment
% ISI 900 ms
%
% April 2015 Anja Roye 
%

clc
close all hidden;
clear all hidden;
% addpath(genpath('C:\gitlab\Anja\experiment'))

%% Get directories & files
DatFolder      = 'C:\DATA\Anja Roye\kennan_oddball\fast\DAT\'; %result data folder 
ExpFolder      = 'C:\DATA\Anja Roye\kennan_oddball\fast\EXP'; %experimtal setup folder 

% Subj		=	input('Enter Subject No: ','s');
Subj = 1;
% Bl          =	input('Enter Block No: ','s');
Bl = 1;
ExpName		=   ['exp_subj' num2str(Subj) 'bl' num2str(Bl) '.mat']; % experimental file  
DatName		=   ['dat_' 'subj' num2str(Subj) 'bl' num2str(Bl)]; %result data file

% RZ6path		=	'C:\DATA\Anja Roye\kennan_oddball\fast\EXP\EegMobitaRZ6.rcx';
RZ6path		= which('EegMobitaRZ6.rcx');

%% Setup TDT 
[ZBus_1, err(1)] = ZBUS(1);
[RZ6_1, err(2)] = RZ6(1,RZ6path);

err
% [RP2_1, err(3), errstr] = RP2(1,RZ6path)

%% Load Experiment
fname = fullfile(ExpFolder,ExpName);
load(fname);

Ntrl	=	size(Par.Grand_mat,1);

%Resp	=	nan(Ntrl,NsigPts);

Q		= struct([]);

for iTrl=1:Ntrl
	if iTrl ==1
		disp('-------------------------');
		%disp('-silence interval 30 sec-');
		%pause(30);
	end
	
	%-- Initialize --%
	NsigPts = settrial(RZ6_1,Par,iTrl);
    
    ZBus_1.zBusTrigA(0,0,6);
	%waits for ISI and checks whether any exit key (spacebar) was pressed
    ch = getkeywait(Par.Grand_mat(iTrl,6)*10^-3); 
    %ch = getkeywait(0.15); 
    if ch == 32
        disp('---------------------------------------------');
		disp(['Trial: ' num2str(iTrl) ' of ' num2str(Ntrl)]);
		disp('ESCAPED');
		cd(DatFolder);
        save([DatName '-' datestr(clock,30)],'Q');
        break
    end
    %alternative
    %pause(Par.Grand_mat(iTrl,6)*10^-3 )
	
	while( RZ6_1.GetTagVal( 'Active') )
		pause(.0001)
	end
	
	%Resp(iTrl,:)	=	RZ6_1.ReadTagV('ADCBuff',0,NsigPts);
	%Trg(iTrl,:)		=	RZ6_1.ReadTagV('TrigBuff',0,NsigPts);
	%Mic(iTrl,:)		=	RZ6_1.ReadTagV('MicBuff',0,NsigPts);
	
	Fs		=	RZ6_1.GetSFreq;
	t		=	linspace(0,NsigPts,NsigPts) / Fs;
	
%	idx		=	find(Resp(iTrl,:) > 0,1,'first');
	
	%figure(1)
	clf

	

    disp('---------------------------------------------');
	disp(['Trial: ' num2str(iTrl) ' of ' num2str(Ntrl)]);
	%disp(['Reaction time: ' num2str(t(idx))]);
    Q(iTrl).subject           = Subj;
    Q(iTrl).condnr			  = Par.Grand_mat(iTrl,3);
	Q(iTrl).freq			  = Par.Grand_mat(iTrl,4);
    Q(iTrl).amp				  = Par.Grand_mat(iTrl,5);
    Q(iTrl).isi				  = Par.Grand_mat(iTrl,6);
%    Q(iTrl).index	 		  =	idx;	%response index 
%	Q(iTrl).reactiontime	  = t(idx); %reaction time in seconds 
	Q(iTrl).fs				  =	Fs;	
	Q(iTrl).t				  = t;
	Q(iTrl).NsigPts			  =	NsigPts;
%	Q(iTrl).response		  = Resp(iTrl,:); 
%	Q(iTrl).trigger			  =	Trg(iTrl,:);	
    
    %% Save data after every trial
	%cd(DatFolder);
	%save([DatName '-' datestr(clock,30)],'Q');
end

if iTrl == Ntrl
	%cd(DatFolder);
	save([DatFolder 'Complete-' DatName '-' datestr(clock,30)],'Q'); %saver after session/completion
end

RZ6_1.Halt;


%RZ6_1.Run;

%function NsigPts = settrial(RP2,PA5_1,PA5_2,Par,TrlNbr)
function NsigPts = settrial(RZ6_1,Par,TrlNbr)
%-- Get info --%
ToneDur		=	Par.ToneDur;
LoTime		=	Par.LoTime;
TriggerDur  =   Par.TriggerDur;
%Npulse		=	Par.Npulse;

Freq		=	Par.Grand_mat(TrlNbr,4);

%Amp = 20; 
%Amp = [35 25 15 5;48 38 28 18;49 39 29 19];

tdtamp = 1; % set intensity

Ear			=	0; %output to both ears, 1..left, 2..right
ISI			=	Par.Grand_mat(TrlNbr,6);
if Freq == 1000
    Npulse=1;
    %TriggerDur=100;
elseif Freq == 1500
    Npulse=2;
    %TriggerDur=50;
end
%Trigger=TrlNbr;
%TriggerDur=0;

%-- Set TDT --%
%RZ6_1.SetTagVal('Trig',Trigger );         
RZ6_1.SetTagVal( 'HiTime', ToneDur  );
RZ6_1.SetTagVal( 'LoTime', LoTime );
RZ6_1.SetTagVal( 'NPulse', Npulse );
RZ6_1.SetTagVal( 'Frequency', Freq );
RZ6_1.SetTagVal( 'Amplitude', tdtamp );
RZ6_1.SetTagVal( 'ISI', ISI );
RZ6_1.SetTagVal( 'Thi', TriggerDur  );

if( Ear == 0 )
	RZ6_1.SetTagVal( 'Left', 1 );
	RZ6_1.SetTagVal( 'Right', 1 );
elseif( Ear == 1 )
	RZ6_1.SetTagVal( 'Left', 1 );
	RZ6_1.SetTagVal( 'Right', 0 );
elseif( Ear == 2 )
	RZ6_1.SetTagVal( 'Left', 0 );
	RZ6_1.SetTagVal( 'Right', 1 );
end

%-- Set the recording buffer --%
Fs			=	RZ6_1.GetSFreq;

NsigPts		=	ceil((ToneDur+LoTime+ISI)*Fs/1000);
RZ6_1.WriteTagV( 'BuffSize', 0, NsigPts );
