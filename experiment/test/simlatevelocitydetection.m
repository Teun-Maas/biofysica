function tmp
close all;
% clearvars;


Fs			= 25000;
dur			= 1;
nsamples	= dur*Fs;
t			= linspace(0,dur,nsamples);
x			= 30*logisticfun(t,0.3,0.05,0.1)+1;
% x			= 30*logisticfun(t,0.3,0.05,0.1)+0; % very noisy velocity profile near 0
% x			= 30*logisticfun(t,0.3,0.05,0.1)-1; % velocity peak when crossing thorugh 0

x = x+0.1*randn(size(x)); % add some noise
x = lowpass(x,'Fc',80,'Fs',Fs,'order',10); % remove some noise
v = [0 diff(x)]*Fs; % determine velocity


subplot(211)
plot(t,x)
nicegraph
xlabel('Time (ms)');
ylabel('Position (deg)');

subplot(212)
plot(t,v)
nicegraph
xlabel('Time (ms)');
ylabel('Velocity (deg/s)');

%%
r = NaN(size(x));
for ii = 1:nsamples
	sel = t>t(ii)-5/1000 & t<t(ii);
	r(ii) = rms(x(sel)).*sign(x(ii)); % remove some noise
end
vr = [0 diff(r)]*Fs; % determine velocity

subplot(211)
hold on
plot(t,r);


subplot(212)
hold on
plot(t,vr);

%%
r = x;
shft = 100;
r1 = circshift(r,1*shft);
r2 = circshift(r,2*shft);
r3 = circshift(r,3*shft);
r4 = circshift(r,4*shft);
r5 = circshift(r,5*shft);

% vdelay = r-r1+r1-r2;
vdelay = (r-r1)/(shft)*Fs;
% vdelay = (r-r5)/(shft*5)*Fs;

vdelay = vdelay;
subplot(211)
hold on
plot(t,r5);

subplot(212)
hold on
plot(t,vdelay,'LineWidth',2);
% vdelay

ylim([-100 1000]);