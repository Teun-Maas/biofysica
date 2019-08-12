fsRZ6               = 48828;

%TRGsel:
%Bit 0 = ZBusB
%Bit 1 = External
%Bit 2 = Software

%SNDsel

%N.B. in case of WAV only two different WAV files (WAV-A for channel A and WAV-B for channel B) can be played during a trial. But they can be played multiple times.
%N.B. in case of A=B the SND-A or WAV-A is played at channel A and channel B

%ACQsel:
%Bit 0..7 = Acquisition Ch0..Ch7



MUX0                = 1*16;
MUX1                = 3*16;
MUX2                = 2*16;
MUX3                = 3*16;
MUX_Set             = 64;
MUX_Rst             = 128;

SPKMOV(1)           = MUX0 + 0;  %mux + index
SPKMOV(2)           = MUX1 + 0;
SPKMOV(3)           = MUX0 + 1;
SPKMOV(4)           = MUX1 + 1;
SPKMOV(5)           = MUX0 + 2;
SPKMOV(6)           = MUX1 + 2;
SPKMOV(7)           = MUX0 + 3;
SPKMOV(8)           = MUX1 + 3;
SPKMOV(9)           = MUX0 + 4;
SPKMOV(10)          = MUX1 + 4;
SPKMOV(11)          = MUX0 + 5;
SPKMOV(12)          = MUX1 + 5;
SPKMOV(13)          = MUX0 + 6;
SPKMOV(14)          = MUX1 + 6;
SPKMOV(15)          = MUX0 + 7;
SPKMOV(16)          = MUX1 + 7;
SPKMOV(17)          = MUX0 + 8;
SPKMOV(18)          = MUX1 + 8;
SPKMOV(19)          = MUX0 + 9;
SPKMOV(20)          = MUX1 + 9;
SPKMOV(21)          = MUX0 + 10;

MOV_Nr_of_Ch        = 21;
MOV_Sp0_is_A        = 1;

start               = 1;
stop                = 0;

no_mov              = 0;
sine_mov            = 1;

A0                  = 1;
A1                  = 2;
A2                  = 4;
A3                  = 8;
A4                  = 16;
A5                  = 32;
A6                  = 64;
A7                  = 128;

B0                  = 1;
B1                  = 2;
B2                  = 4;
B3                  = 8;
B4                  = 16;
B5                  = 32;
B6                  = 64;
B7                  = 128;

STM_nTrials         = 6;
STM_nPars           = 7;

LPfreq              = 10000;
HPfreq              = 300;

Centre_Freq         = 500;       
Mod_Freq            = 1;         
Mod_BW              = 25;         

TRG_ZBusB           = 1;
INP_EVT             = 2;
TRG_Soft            = 4;

TRG_Mask            = TRG_ZBusB +INP_EVT;   
EVT_Mask            = A4;
INP_Mask            = A0 + A3;

TSK_Wait_For_TRG    = 0;
TSK_SND_A           = 1;
TSK_SND_B           = 2;
TSK_MUX             = 3;
TSK_Signaling       = 4;
TSK_MOV             = 5;
TSK_ACQ             = 6;
TSK_DIO_out         = 7;
TSK_TRG_Out         = 8;
TSK_Reset           = 9;
TSK_Ready           = 10;
TSK_Hold_INP        = 11;

SND_None            = 0;
SND_Stop            = 1;
SND_Tone            = 2;
SND_Sweep           = 3;
SND_Noise           = 4;
SND_Ripple          = 5;
SND_WAV             = 6;
SND_B_is_A          = 7;



% construction of STM data matrix
STM_1(1)            = TSK_Wait_For_TRG; % TSK type
STM_2(1)            = 0;         % sound type none
STM_3(1)            = 0;         %delay in milisecs
STM_4(1)            = TRG_Mask;   %Par1
STM_5(1)            = EVT_Mask;        %Par2
STM_6(1)            = 0;         %Par3
STM_7(1)            = 0;         %Par4 %attenuation A&B

STM_1(2)            = TSK_SND_A;        % TSK type
STM_2(2)            = SND_Tone;  % sound type
STM_3(2)            = 1000;      %delay in milisecs
STM_4(2)            = Centre_Freq;       %Par1
STM_5(2)            = Mod_Freq;         %Par2 
STM_6(2)            = Mod_BW;         %Par3
STM_7(2)            = 5;        %Par4 %attenuation A&B

STM_1(3)            = TSK_Hold_INP;        % TSKion type
STM_2(3)            = 0;  % sound type
STM_3(3)            = 1000;      %delay in milisecs
STM_4(3)            = 1;%Par1
STM_5(3)            = INP_Mask;         %Par2
STM_6(3)            = 0;         %Par3  ACQ Ch1
STM_7(3)            = 0;        %Par4 %attenuation A&B

STM_1(4)            = TSK_Wait_For_TRG;        % TSKion type
STM_2(4)            = 0;
STM_3(4)            = 2001;      % delay
STM_4(4)            = TRG_Mask;    % ACQ Ch1  ACQ start
STM_5(4)            = EVT_Mask;       %Par2 ACQ time
STM_6(4)            = 0;     %Par3 ACQ Ch1
STM_7(4)            = 0;        %Par4 %attenuation A&B

STM_1(5)            = TSK_SND_A;        % TSKion type
STM_2(5)            = SND_Stop;      % sound type
STM_3(5)            = 0;      %delay in milisecs
STM_4(5)            = 0;         %Par1
STM_5(5)            = 0;         %Par2
STM_6(5)            = 0;         %Par3
STM_7(5)            = 0;        %Par4 %attenuation A&B

STM_1(6)            = TSK_Hold_INP;% TSKion type
STM_2(6)            = 0;         % sound type
STM_3(6)            = 1000;      %delay in milisecs
STM_4(6)            = 0;         %Par1
STM_5(6)            = 0;         %Par2
STM_6(6)            = 0;         %Par3
STM_7(6)            = 0;        %Par4 %attenuation A&B


STM_Matrix = [STM_1(1) STM_2(1) STM_3(1) STM_4(1) STM_5(1) STM_6(1) STM_7(1),
              STM_1(2) STM_2(2) STM_3(2) STM_4(2) STM_5(2) STM_6(2) STM_7(2),
              STM_1(3) STM_2(3) STM_3(3) STM_4(3) STM_5(3) STM_6(3) STM_7(3),
              STM_1(4) STM_2(4) STM_3(4) STM_4(4) STM_5(4) STM_6(4) STM_7(4),
              STM_1(5) STM_2(5) STM_3(5) STM_4(5) STM_5(5) STM_6(5) STM_7(5),
              STM_1(6) STM_2(6) STM_3(6) STM_4(6) STM_5(6) STM_6(6) STM_7(6)];


% communication with RZ6

zBus = actxcontrol('ZBUS.x',[1 1 1 1]);

if zBus.ConnectZBUS('GB')
 e = 'ZBus connected';
else
 e = 'Unable to connect ZBus'
end  

%Connects to RZ6 #1 via Optical Gigabit
RZ6=actxcontrol('RPco.x',[5 5 26 26]);

if RZ6.ConnectRZ6('GB',1)
 e = 'RZ6 connected';
else
 e = 'Unable to connect RZ6'
end

if RZ6.LoadCOF('RCX_Uni_1C_50kHz_V3.00.rcx')
 e='LoadCOF OK';
else
 e='LoadCOF error'
end

RZ6.Run;

if RZ6.WriteTagVEX('STM_Matrix', 0, 'I32', transpose(STM_Matrix)) %transpose to single column
 e='WriteTagVex OK';
else
 e='WriteTagVEX1 error'
end

if RZ6.WriteTagVEX('MOV_Sp_Array', 0, 'I32', transpose(SPKMOV)) %transpose to single column
 e='WriteTagV OK';
else
 e='WriteTagVEX2 error'
end

if RZ6.SetTagVal('MOV_Sp0_is_A', MOV_Sp0_is_A ) %transpose to single column
 e='SetTagV OK';
else
 e='SetTagVal error'
end


pause(0.1);

Matrix1 = RZ6.ReadTagVEX('STM_Matrix',0,STM_nTrials,'I32','I32',STM_nPars)
%Matrix2 = RZ6.ReadTagVEX('MOV_Sp_Array',0, 1,'I32','I32',MOV_Nr_of_Ch)




