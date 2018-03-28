close all
clearvars


%% Unaided
x = -90:90;
y = 0.4*x+50+30*randn(size(x));
y(y>70) = 70;

sel = x<0;
MAE(1) = mean(abs(y(sel)-x(sel)));
sel = x>0;
MAE(2) = mean(abs(y(sel)-x(sel)));

subplot(121)
plot(x,y,'ko','MarkerFaceColor','w');

xlim([-90 90]);
ylim([-90 90]);
nicegraph

title(MAE)



%% Aided
x = -90:90;
y = 0.6*x+20+30*randn(size(x));
y(y>70) = 70;

sel = x<0;
MAE(1) = mean(abs(y(sel)-x(sel)));
sel = x>0;
MAE(2) = mean(abs(y(sel)-x(sel)));

subplot(122)
plot(x,y,'ko','MarkerFaceColor','w');

xlim([-90 90]);
ylim([-90 90]);
nicegraph

title(MAE)
