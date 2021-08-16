
times=0:0.1:60;
y1=merfeld(90, 10, times, 5, 5);
y2=zeros(size(y1));

plot(times, y1, '.');
figure(gcf);

vs=vs_servo;

vs.write_profile(y1,y2);
delete(vs);



