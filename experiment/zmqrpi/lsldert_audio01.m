client=lsldert_client('raspi4');

fsamp=44100;
fsine=1000;
duration=2;

nsamp=fsamp*duration;
times=(0:nsamp-1)*duration/fsamp;
pcm_data=int16(32767*sin(2*pi*fsine*times));
pcm_header=uint16([0, 0, 0, 0, 0, 0]);
pcm_left=pcm_data';
[~,telapsed]=client.send('S FILL 1', pcm_header, pcm_data);
telapsed
delete(client);
