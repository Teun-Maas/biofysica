function genexp_defaultcal
% GENEXP_DEFAULTCAL
%
% This will generate an EXP-file for a default calibration experiment. 
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

expfile     = 'defaultcal.exp'; % file name
datdir      = 'DEFAULT';
%% Some Flags
showexp     = true;

%% Desired azimuth and elevation
% des_az		= -90:15:90;
% des_el		= -90:15:90;
des_az		= -90:30:90;
des_el		= [-45 -90:30:90];
[des_az,des_el] = meshgrid(des_az,des_el);
des_az		= des_az(:);
des_el		= des_el(:);
sel			= (abs(des_az)+abs(des_el))<=90 & des_el>-60; 
des_az		= des_az(sel);
des_el		= des_el(sel);

%% Actual azimuth and elevation
% The actual speaker positions are not perfectly aligned with 15 deg
% Let's see what they are
cfg			= spherelookup; % sphere positions
channel		= cfg.interpolant(des_az',des_el');
X		= cfg.lookup(channel+1,5);
Y		= cfg.lookup(channel+1,6);

%% Present systematically (from left to right, from bottom to up)
XY = [des_az des_el];
[~,idx] = sortrows(XY,[1 2]);
X = X(idx);
Y = Y(idx);

%% Graphics

close all
plot(des_az,des_el,'.')
hold on
plot(X,Y,'.')

axis([-130 130 -130 130]);
axis square
set(gca,'TickDir','out');
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');


%% Number and size
Sz				= size(X); 
N				= Sz(1);% number of trials


%% Save data somewhere
writeexp(expfile,datdir,X,Y); 
% see below, these are helper functions to write an exp-file line by line / stimulus by stimulus

%% Show the exp-file in Wordpad
% for PCs
if ispc && showexp
    dos(['"C:\Program Files\Windows NT\Accessories\wordpad.exe" ' expfile ' &']);
end

function writeexp(expfile,datdir,theta,phi)
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
	writeled(fid,'LED',round(theta(ii)),round(phi(ii)),5,0,0,1,100); % fixation LED for 100 ms after button press
	writetrg(fid,1,2,0,0,1);		% Button trigger after LED has been fixated
	writeacq(fid,1,0);	% Data Acquisition immediately after button press
	
end
fclose(fid);

