t=0:0.01:20;

y=sinetrpz(45, 2, t, 2, 5);
plot(t,y,'.');
figure(gcf);