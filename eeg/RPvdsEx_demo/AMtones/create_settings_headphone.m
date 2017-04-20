NTrials = 24;
Saving_dir  = 'C:\Users\Elisabeth Noordanus\Documents\MATLAB\';  % End with backslash. Dir in which to save the produced file 
% Matlab always has write access to the MATLAB folder in C:\Users\ username \Documents\MATLAB\
% If you want to look at/change the file: double click it in Matlab,
% Windows has a problem with recognizing the correct file type if you
% didn't change the file association for the .mat extension.
Subj    	= input('Enter the subject code for the data file: ','s');
Bl      	= input('Enter the block number for the data file: ','s');

%-- Intialization of variables that don't change during the experiment --%
Att_A           = 60;       % Attenuation (db) of the left output
Att_B           = 60;       % Attenuation (db) of the right output
Stim_time_l     = 1300;     % Duration in msec when only the left sound is played
Stim_time_r     = 1300;     % Duration in msec when only the right sound is played
Stim_time_both  = 12300;    % Duration in msec when both sounds are played
A_m_l           = 2.5;      % Amplitude (V) of the modulation of the left tone
A_m_r           = 2.5;      % Amplitude (V) of the modulation of the right tone
Ashift_m_l      = 2.5;      % Shift of the amplitude (V) of the modulation of the left tone
Ashift_m_r      = 2.5;      % Shift of the amplitude (V) of the modulation of the right tone
Ashift_c_l      = 0;        % Shift of the amplitude (V) of the carrier wave of the left tone
Ashift_c_r      = 0;        % Shift of the amplitude (V) of the carrier wave of the right tone
f_m_l           = 37;       % Frequency (hz) of the modulation of the left tone
f_m_r           = 43;       % Frequency (hz) of the modulation of the right tone
f_c_l           = 1000;     % Frequency (hz) of the carrier wave of the left tone
f_c_r           = 2000;     % Frequency (hz) of the carrier wave of the right tone
phase_m_l       = 0;        % Starting phase (rad) of the modulation of the left tone
phase_m_r       = 0;        % Starting phase (rad) of the modulation of the right tone
phase_c_l       = 0;        % Starting phase (rad) of the carrier wave of the left tone
phase_c_r       = 0;        % Starting phase (rad) of the carrier wave of the right tone

% The settings to test
% Take care that NTrials / # of variables gives an integer value, otherwise no
% equal amount of presentations
var_array1 = [1,2]; % left = 1 / right = 2.
settings = [];
N=size(var_array1,2);
for teller = 1:N
    settings = [settings; (var_array1(teller) * ones(round(NTrials/N),1))];
end
if (size(settings) > NTrials)  % In case NTrials / # of variables is no integer
    settings = settings(1:NTrials,:);
elseif (size(settings,1) > NTrials)
    while (size(settings,1) < NTrials)
        settings = [settings; settings(end,:)];
    end
end
% array presented_stimuli consists of: [ISI_1,ISI_2,left / right]
% ISI_1: Time (sec) to wait before presenting the sound which indicates on which side to concentrate
% ISI_2: Time (sec) to wait between presenting the sound which                             
presented_stimuli_all=[(4.5+rand(NTrials,1)), ones(NTrials,1), settings];
% Shuffle all trials
idx = randperm(NTrials);
presented_stimuli = presented_stimuli_all(idx,:);
% Make the structure with all the settings
exp_settings = struct('Att_A',Att_A,'Att_B', Att_B,'Stim_time_l', Stim_time_l,'Stim_time_r', Stim_time_r, ...
    'Stim_time_both', Stim_time_both,'A_m_l', A_m_l,'A_m_r', A_m_r,'Ashift_m_l', Ashift_m_l,'Ashift_m_r', Ashift_m_r, ...
    'Ashift_c_l', Ashift_c_l,'Ashift_c_r', Ashift_c_r,'f_m_l', f_m_l,'f_m_r', f_m_r,'f_c_l', f_c_l, ...
    'f_c_r', f_c_r,'phase_m_l', phase_m_l,'phase_m_r', phase_m_r,'phase_c_l', phase_c_l,'phase_c_r', phase_c_r, ...
    'presented_stimuli', presented_stimuli);

% Save to file
FileName=[Saving_dir, 'AMtones-headphone_subj_',Subj,'_block_',Bl,'.mat'];
save(FileName, 'exp_settings');

%{
% Make test version of the file, write correct load dir and filename below.
% Fill in also correct name of the saved file
% Write in the idx line the correct # of trials. In the test block 8 trials
% are saved.
% Check and correct if not enough semi random.
load('C:\Matlab testen RU\EEG\EEG exp 1\expblock\AMtones-headphone_subj_1_block_3.mat');
length = size(exp_settings.presented_stimuli,1);
idx = randperm(length);
exp_settings.presented_stimuli = exp_settings.presented_stimuli(idx,:);
exp_settings.presented_stimuli = exp_settings.presented_stimuli(1:8,:);
save('C:\Users\Elisabeth Noordanus\Documents\MATLAB\AMtones-headphone_subj_test1_block_test1.mat', 'exp_settings');
%}

%{
% Make speakers version from headphone mat
length = size(exp_settings.presented_stimuli,1);
exp_settings.Att_A = 0;
exp_settings.Att_B = 0;
exp_settings.presented_stimuli(:,4)=exp_settings.presented_stimuli(:,3);
exp_settings.presented_stimuli(:,3) = exp_settings.presented_stimuli(:,1)-2;
exp_settings.presented_stimuli(1:length/2,1) = 37;
exp_settings.presented_stimuli(1:length/2,2) = 1000;
exp_settings.presented_stimuli(length/2+1:end,1) = 43;
exp_settings.presented_stimuli(length/2+1:end,2) = 2000;
idx = randperm(length);
exp_settings.presented_stimuli = exp_settings.presented_stimuli(idx,:);
save('C:\Users\Elisabeth Noordanus\Documents\MATLAB\AMtone-speaker_subj_test1_block_test1-nieuw.mat', 'exp_settings');

% Change general settings
exp_settings.Att_A = 60;
exp_settings.Att_B = 60;
exp_settings.A_m_l = 2.5;
exp_settings.A_m_r = 2.5;
exp_settings.Ashift_m_l = 2.5;
exp_settings.Ashift_m_r = 2.5;

%}