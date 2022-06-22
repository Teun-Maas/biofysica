function genexp_defaultloc
% GENEXP_DEFAULTLOC
%
% This will generate an EXP-file for a default localization experiment. 
% EXP-files are used for the psychophysical experiments at the
% Biophysics Department of the Donders Institute for Brain, Cognition and
% Behavior of the Radboud University Nijmegen, the Netherlands.
% 
% See also the manual for the experimental set-ups at www.mbys.ru.nl/staff/m.vanwanrooij/neuralcode/manuals/audtoolbox.pdf.
% See also WRITESND, WRITELED, WRITETRG, GENWAV_DEFAULT, etc


%% Initialization
home;
close all;
clear hidden;
disp('>> GENERATING EXPERIMENT <<');

%% Default input
% Minimum and maximum fixation light offset after button press
minled      = 300;
maxled      = 800;
expfile     = 'defaultloc.exp'; % file name
datdir      = 'DEFAULT';
% Minimum and maximum sound onset after LED offset
minsnd      = 200; % 200 ms gap (total darkness and silence) for optimal OPN inhibition...
maxsnd      = 200;
%% Some Flags
showexp     = true;

%% Desired azimuth and elevation
des_az		= -90:15:90;
des_el		= -90:15:90;
% des_az		= -90:30:90;
% des_el		= [-50 -30:30:90];
[des_az,des_el] = meshgrid(des_az,des_el);
des_az		= des_az(:);
des_el		= des_el(:);
sel			= (abs(des_az)+abs(des_el))<=90 & des_el>-60; 
des_az		= des_az(sel);
des_el		= des_el(sel);
nloc		= numel(des_az);

% choose nloc random locations
% n			= numel(des_az);
% nloc		= 30;
% idx			= randperm(n,nloc);
% des_az		= des_az(idx);
% des_el		= des_el(idx);

%% Actual azimuth and elevation
% The actual speaker positions are not perfectly aligned with 15 deg
% Let's see what they are
cfg			= spherelookup; % sphere positions
channel		= round(cfg.interpolant(des_az',des_el'));
% channel+1
X		= cfg.lookup(channel+1,5);
Y		= cfg.lookup(channel+1,6);

%% Graphics
% az = cfg.lookup(:,5);
% el = cfg.lookup(:,6);

close all
plot(des_az,des_el,'.')
hold on
plot(X,Y,'.')

axis([-130 130 -130 130]);
axis square
set(gca,'TickDir','out');
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');

%% Intensities & frequency bands
int				= [45 55 65]; % approx. dB
freq			= [1 2 3]; % BB, HP, LP
[X,~,~]			= ndgrid(X,int,freq);
[Y,int,freq]	= ndgrid(Y,int,freq);
X				= X(:);
Y				= Y(:);
int				= int(:);
freq			= freq(:);

%% Number and size
Sz				= size(X); 
N				= Sz(1);% number of trials


%% Randomize sound samples (to simulate fresh noise
snd				= freq;
sel				= snd==1; % BB
p				= transpose(randperm(100,sum(sel))-1);
snd(sel)		= snd(sel)*100 + p;

sel				= snd==2; % HP
p				= transpose(randperm(100,sum(sel))-1);
snd(sel)		= snd(sel)*100 + p;

sel				= snd==3; % LP
p				= transpose(randperm(100,sum(sel))-1);
snd(sel)		= snd(sel)*100 + p;

%% Get random timings
% Again, randomization is a good thing. When does the LED go off after
% button press, when does the sound go off after LED extinction?
ledon			= rndval(minled,maxled,Sz);
% Choose a value that hinders the subjects to enter in a default mode.
sndon			= rndval(minsnd,maxsnd,Sz); 
% By default 200 ms, this gap is usually sufficient to reduce fixation inhibition of reaction times of saccades. 

%% Randomize
rnd				= randperm(N);
X				= round(X(rnd));
Y				= round(Y(rnd));
ledon			= ledon(rnd);
sndon			= sndon(rnd);
int				= int(rnd);
snd				= snd(rnd);


%% Save data somewhere
writeexp(expfile,datdir,X,Y,snd,int,ledon,sndon); 
% see below, these are helper functions to write an exp-file line by line / stimulus by stimulus

%% Show the exp-file in Wordpad
% for PCs
if ispc && showexp
    dos(['"C:\Program Files\Windows NT\Accessories\wordpad.exe" ' expfile ' &']);
end

function writeexp(expfile,datdir,theta,phi,snd,int,ledon,sndon)
% Save known trial-configurations in exp-file
%
%WRITEEXP WRITEEXP(FNAME,DATDIR,THETA,PHI,ID,INT,LEDON,SNDON)
%
% WRITEEXP(FNAME,THETA,PHI,ID,INT,LEDON,SNDON)
%
% Write exp-file with file-name FNAME.
%
%
% See also manual at neural-code.com
expfile		= fcheckext(expfile,'.exp'); % check whether the extension exp is included


fid         = fopen(expfile,'w'); % this is the way to write date to a new file
ntrials     = numel(theta); % only 135 trials

%% Header of exp-file
ITI			= [0 0];  % useless, but required in header
Rep			= 1; % we have 0 repetitions, so insert 1...
Rnd			= 0; % we randomized ourselves already
Mtr			= 'n'; % the motor should be on
writeheader(fid,datdir,ITI,ntrials*Rep,Rep,Rnd,Mtr,'Lab',2); % helper-function

%% Body of exp-file
% Create a trial
for ii               = 1:ntrials		% each location
	writetrl(fid,ii);
	writeled(fid,'LED',0,0,5,0,0,1,ledon(ii)); % fixation LED
	
	writetrg(fid,1,2,0,0,1);		% Button trigger after LED has been fixated
	writeacq(fid,1,ledon(ii));	% Data Acquisition immediately after fixation LED exinction
	
	writesnd(fid,'SND',round(theta(ii)),phi(ii),snd(ii),int(ii),1,ledon(ii)+sndon(ii)); % Sound on
end
fclose(fid);

