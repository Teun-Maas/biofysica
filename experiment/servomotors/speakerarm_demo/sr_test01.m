
times=0:0.1:60;
sv_vert=sinetukey(100, 10, times, 5);
sv_chair=zeros(size(sv_vert));
sv_horiz=zeros(size(sv_chair));
sv_horiz=sv_vert;

plot(times, sv_vert, '.');
figure(1);

sr=sr_servo;

sr.write_profile(sv_vert, sv_chair, sv_horiz);
sr.enable();
input('press enter to start')
sr.start();
input('press enter to stop')
sr.stop();
[pv_vert, pv_chair, pv_horiz]=sr.read_profile_pv();
figure(2);
plot(times,pv_vert(1:size(times,2)));
delete(sr);
