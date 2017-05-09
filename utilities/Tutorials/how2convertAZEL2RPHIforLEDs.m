close all
clearvars;


phi = 0:360;
r = repmat(10,size(phi));
[azel] = pa_rphi2azel(r,phi);
plot(azel(:,1),azel(:,2));
hold on
r = repmat(30,size(phi));
[azel] = pa_rphi2azel(r,phi);
plot(azel(:,1),azel(:,2));
r = repmat(45,size(phi));
[azel] = pa_rphi2azel(r,phi);
plot(azel(:,1),azel(:,2));

r = repmat(60,size(phi));
[azel] = pa_rphi2azel(r,phi);
plot(azel(:,1),azel(:,2));
axis square
axis([-90 90 -90 90])
% horline(-75:15:75);
% verline(-75:15:75);
% horline(-75:15:75);
% verline(-75:15:75);
set(gca,'XTick',[-90:15:-30 -15:5:15 30:15:90]);
set(gca,'YTick',[-90:15:-30 -15:5:15 30:15:90]);
grid on

az = [-90:15:-30 -15:5:15 30:15:90];
[az,el] = meshgrid(az,az);
az = az(:);
el = el(:);
sel = (abs(az)+abs(el))<=90;
az = az(sel);
el = el(sel);
[r,phi] = pa_azel2pol(az,el);
sel = r<=45;
az = az(sel);
el = el(sel);

plot(az,el,'o');

title(sum(sel));