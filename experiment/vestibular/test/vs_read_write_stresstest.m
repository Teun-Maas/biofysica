vs=vs_servo;

n=2000;

for ii=1:100
    times=(0:n-1)/10;
    svax1 = zeros(1,n);
    svax2 = zeros(1,n);
    fprintf('write: ');
    tic;
    vs.write_profile(svax1,svax2);
    toc
    tic;
    fprintf('read: ');
    [svax1, svax2]=vs.read_profile_sv;
    [pvax1, pvax2]=vs.read_profile_pv;
    toc
end
delete(vs);


clf;
plot(times,svax2,'-');
hold('on');
%plot(times,pvax2,'.');
grid('on');
xlabel('time [s]');
ylabel('angle [deg]');
figure(gcf);

