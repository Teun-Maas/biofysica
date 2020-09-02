client=lsldert_pubclient('raspi6');

% fsamp=44100;
% fsine=1000;
% duration=2;
% nsamp=fsamp*duration;
% times=(0:nsamp-1)*duration/fsamp;
% pcm_data=single(sin(2*pi*fsine*times));

% these are example sounds from matlab, mono, 8192 S/s, Nx1 matrix
%load handel
load laughter
%load gong
%load train
%load chirp
%load splat

[y, Fs]=audioread('mariobiondi.mp3',[10 20]*44100);
% let it be mono...
pcm_data=y;
% ...or make it two channel N by 2 matrix
%pcm_data=[y,y];
% --
nsamp=max(size(pcm_data));
nchan=min(size(pcm_data));
pcm_header=uint32([Fs, nchan, nsamp, 0, 0, 0]);
% y2=[y';y'];
% pcm_data=reshape(y2,1,[]);
% nchan=2;
% nsamp=0.5*max(size(pcm_data));
% pcm_header=uint32([Fs, nchan, nsamp, 0, 0, 0]);
% %pcm_data(:)=0;
[~,telapsed]=client.send('AF 0', pcm_header, single(pcm_data));
[~,telapsed]=client.send('AP 0');
telapsed
input('press enter >>>>');
%[~,telapsed]=client.send('AS 0');

delete(client);
