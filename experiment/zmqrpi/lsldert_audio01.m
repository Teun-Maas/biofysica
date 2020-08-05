client=lsldert_pubclient('raspi6');

% fsamp=44100;
% fsine=1000;
% duration=2;
% nsamp=fsamp*duration;
% times=(0:nsamp-1)*duration/fsamp;
% pcm_data=single(sin(2*pi*fsine*times));

load handel
y2=[y',y'];
y=reshape(y2,[],1);
nchan=2;
nsamp=0.5*max(size(y));
pcm_data=y;
pcm_header=uint32([Fs, nchan, nsamp, 0, 0, 0]);

pcm_left=pcm_data';
[~,telapsed]=client.send('AF 3', pcm_header, single(pcm_data));
[~,telapsed]=client.send('AP 3');
telapsed
delete(client);
