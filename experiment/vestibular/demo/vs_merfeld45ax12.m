
times=0:0.1:60;
y1=merfeld(90, 10, times, 5, 5);
y2=merfeld(30, 15, times, 5, 5);

plot(times, y1,times,y2, '.');
figure(gcf);

vs=vs_servo;

vs.write_profile(y1,y2);
delete(vs);



