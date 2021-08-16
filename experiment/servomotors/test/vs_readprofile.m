vs=vs_servo;


times=(0:1999)/10;
[svax1, svax2]=vs.read_profile_sv;
[pvax1, pvax2]=vs.read_profile_pv;
delete(vs);


clf;
plot(times,svax2,'-');
hold('on');
plot(times,pvax2,'.');
grid('on');
xlabel('time [s]');
ylabel('angle [deg]');
figure(gcf);

