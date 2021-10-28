
times=0:0.1:20;
vert=sinetukey(100, 10, times, 1);
%chair=zeros(size(vert));
%horiz=zeros(size(chair));
horiz=vert;
chair=vert;


sr=sr_servo;

sr.write_profile(vert, chair, horiz);
sr.enable();
input('press enter to start')
sr.start();
input('press enter to stop')
sr.stop();
[pv_vert, pv_chair, pv_horiz]=sr.read_profile_pv();
[sv_vert, sv_chair, sv_horiz]=sr.read_profile_sv();

range=1:length(times);

%% Vertical
figure('Name','Vert Axis');
title('SV');
subplot(131)
plot(times, vert(range),'.');
title('SV');

subplot(132);
plot(times,pv_vert(range),'.');
title('PV PLC');

subplot(133);
plot(times,sv_vert(range),'.');
title('SV PLC');

%% Horizontal
figure('Name','Horiz Axis');
title('SV');
subplot(131)
plot(times, horiz(range),'.');
title('SV');

subplot(132);
plot(times,pv_horiz(range),'.');
title('PV PLC');

subplot(133);
plot(times,sv_horiz(range),'.');
title('SV PLC');

% Seat
figure('Name','Seat Axis');
title('SV');
subplot(131)
plot(times, chair(range),'.');
title('SV');

subplot(132);
plot(times,pv_chair(range),'.');
title('PV PLC');

subplot(133);
plot(times,sv_chair(range),'.');
title('SV PLC');
delete(sr);
