client=lsldert_pubclient('raspi6');

% fsamp=44100;
% fsine=1000;
% duration=2;
% nsamp=fsamp*duration;
% times=(0:nsamp-1)*duration/fsamp;
% pcm_data=single(sin(2*pi*fsine*times));

load handel
y2=[y';y'];
pcm_data=reshape(y2,1,[]);
nchan=2;
nsamp=0.5*max(size(pcm_data));
pcm_header=uint32([Fs, nchan, nsamp, 0, 0, 0]);
pcm_data(:)=0;
[~,telapsed]=client.send('AF 3', pcm_header, single(pcm_data));
[~,telapsed]=client.send('AP 3');
telapsed
pause(4);
[~,telapsed]=client.send('AS 3');

delete(client);
