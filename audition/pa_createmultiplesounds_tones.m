function pa_createmultiplesounds_tones

%%
close all hidden
clear all hidden %#ok<CLFUN>

pa_genexp_tone

freq = pa_oct2bw(1000,-1:4); % Hz
dur = [100 2000]; % ms
[F,D] = ndgrid(freq,dur);
F = F(:);
D = D(:);
N = numel(F)
T = 29*2*9*N*2/60/60 % hours
for indx = 1:N
	
	close all
	
	%% Default
	Fs		= 48828.125; % Freq (Hz)
	Nramp	= round(5/1000*Fs); % 5 ms
	Nlevel	= round(20/1000*Fs); % 5 ms
	eqfile = 'E:\MATLAB\PANDA\Donders\Setups\supplementary\visaton_equalizer.mat';
	
	col = bone(4);
	
	%% Step 1. Generate sound
	% pa_gensweep - Schroeder sweep
	% pa_gengwn - Gaussian white noise from random time trace
	% pa_gengwnflat - Gaussian white noise from random phase spectrum
	% pa_gentone - Pure tone
	% pa_genripple - generate ripple
	
	snd = pa_gentone(F(indx),D(indx));
	snd = snd./max(abs(snd));
	
	subplot(211)
	plot(snd,'k-','Color',col(1,:),'LineWidth',4)
	hold on
	xlim([1 numel(snd)]);
	xlabel('Time (samples)');
	ylabel('Amplitude (au)');
	box off
	title('Entire waveform');
	
	subplot(245)
	plot(snd(1:Nramp*3),'k-')
	hold on
	xlim([1 Nramp*3]);
	xlabel('Time (samples)');
	ylabel('Amplitude (au)');
	title('Step 1. Signal generation');
	box off
% 	p		= audioplayer(snd,Fs);
	% playblocking(p);
	
	
	%% Step 2. Filter
	% pa_highpass
	% pa_lowpass
	% pa_equalizer
	snd = pa_equalizer(snd,eqfile);
	snd = snd./max(abs(snd));
	% snd = pa_lowpass(snd);
	% snd = pa_highpass(snd);
	snd = pa_bandpass(snd);
	subplot(211)
	plot(snd,'k-','Color',col(2,:),'LineWidth',2)
	
	subplot(246)
	plot(snd(1:Nramp*3),'k-')
	xlim([1 Nramp*3]);
	xlabel('Time (samples)');
	ylabel('Amplitude (au)');
	title('Step 2. Equalizer');
	box off
% 	p		= audioplayer(snd,Fs);
% 	playblocking(p);
	
	%% Step 3. Ramp
	% pa_ramp
	snd = pa_ramp(snd,Nramp);
	snd = snd./max(abs(snd));
	
	subplot(211)
	plot(snd,'k-','Color',col(3,:),'LineWidth',.5)
	
	subplot(247)
	plot(snd(1:Nramp*3),'k-')
	xlim([1 Nramp*3]);
	xlabel('Time (samples)');
	ylabel('Amplitude (au)');
	title('Step 3. Ramp');
	box off
% 	p		= audioplayer(snd,Fs);
% 	playblocking(p);
	
	%% Step 4. 20 ms level ramp
	% pa_levelramp
	
	snd = pa_fart_levelramp(snd,Nlevel);
	

	snd = snd./max(abs(snd));
	
	
	subplot(211)
	plot(snd,'k-','Color',col(4,:),'LineWidth',.1)
	
	subplot(248)
	plot(snd(1:Nramp*3),'k-','LineWidth',1)
	xlim([1 Nramp*3]);
	xlabel('Time (samples)');
	ylabel('Amplitude (au)');
	title('Step 4. Zero-prepend');
	box off
% 	p		= audioplayer(snd,Fs);
	% playblocking(p);


	%% Step 5. save
	pa_datadir;
	pa_writewav(snd,['snd' num2str(indx,'%03d')]); % bb
	% pa_writewav(snd,'snd002'); % lp
	% pa_writewav(snd,'snd003'); % hp
	
	drawnow
% 	pause
end
% print('-depsc','-painter',mfilename);


%% 

function pa_genexp_tone
% PA_GENEXP_SPATIALPRIOR
%
% This will generate an EXP-file of a spatial prior learning
% experiment. EXP-files are used for the psychophysical experiments at the
% Biophysics Department of the Donders Institute for Brain, Cognition and
% Behavior of the Radboud University Nijmegen, the Netherlands.
% 
% See also the manual for the experimental set-ups at www.mbys.ru.nl/staff/m.vanwanrooij/neuralcode/manuals/audtoolbox.pdf.
% See also WRITESND, WRITELED, WRITETRG, etc

% (c) 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
home;
disp('>> GENERATING EXPERIMENT <<');
close all
clear all;

% Default input
minled      = 300;
maxled      = 800;
expfile     = 'tonetest0.exp';
datdir      = 'DAT';
minsnd      = 200;
maxsnd      = 200;

%% Some Flags
showexp     = true;

%% Stimulus locations
el = 12;
az = 0;


% transform to hoop coordinates
[theta,phi] = pa_azel2fart(az,el); % theta = hoop azimuth, phi = speaker number, from 1-29 and 101-129, see also manual
theta		= round(theta);

%% By default, sounds are presented at various intensities
int			= 10:10:80; % this is approximately 40, 50 and 60 dBA, 3 levels

%% Now define the IDs of the sounds we want to use.
snd = 1:12;

[snd,int] = meshgrid(snd,int);
snd = snd(:);
int = int(:);


%% Save data somewhere
pa_datadir; % Change this to your default directory
writeexp(expfile,datdir,snd,int); 
% see below, these are helper functions to write an exp-file line by line / stimulus by stimulus

%% Show the exp-file in Wordpad
% for PCs
if ispc && showexp
    dos(['"C:\Program Files\Windows NT\Accessories\wordpad.exe" ' expfile ' &']);
end

function writeexp(expfile,datdir,snd,int)
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
expfile		= pa_fcheckext(expfile,'.exp'); % check whether the extension exp is included


fid         = fopen(expfile,'w'); % this is the way to write date to a new file
ntrials     = numel(snd); % only 135 trials

%% Header of exp-file
ITI			= [0 0];  % useless, but required in header
Rep			= 1; % we have 0 repetitions, so insert 1...
Rnd			= 0; % we randomized ourselves already
Mtr			= 'n'; % the motor should be on
writeheader(fid,datdir,ITI,ntrials*Rep,Rep,Rnd,Mtr); % helper-function
unique(snd)

%% Body of exp-file
% Create a trial
for ii               = 1:ntrials		% each location
	writetrl(fid,ii);
	
	
	writesnd(fid,'SND1',0,12,snd(ii),int(ii),0,100); % Sound on
	writeinp(fid,1);
end
fclose(fid);


