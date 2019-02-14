
close all;
clear all;

HSfun = @(x) 9.7*sin(0.02*x+0.27);
spk1 = 20;
spk2 = 25;
dx = abs(spk2-spk1);
desired_angle	= spk1:0.1:spk2;
desired_angle	= desired_angle(1:end-1);
HS				= HSfun(desired_angle);
figure(1)
clf
subplot(231)
plot(desired_angle,HS);

%establish loudspeaker locations and numbers


%determine panning levels
s1_level = sin(abs(rem(desired_angle,dx)/dx*(pi/2)));
s2_level = cos(abs(rem(desired_angle,dx)/dx*(pi/2)));

% s1_level = s1_level
% s1_level = linspace(0,1,size(s1_level,2));
% s2_level = fliplr(s1_level);
whos s2_level s1_level
% %swap if s1 is lower in level
% [s1_level s2_level] = deal(max([s1_level s2_level]),min([s1_level s2_level]));


subplot(234)
plot(desired_angle,s1_level);
hold on
plot(desired_angle,s2_level);

HSspk1 = HS(1)+60;
HSspk2 = HS(end)+60;

amp1 = (10^(HSspk1/20)*s1_level).^2+(10^(HSspk2/20)*s2_level).^2;
amp2 = (10^(HSspk1/20)*s2_level).^2+(10^(HSspk2/20)*s1_level).^2;
dB1 = 10*log10(amp1);
dB2 = 10*log10(amp2);
ILD = dB2-dB1;

subplot(235)
plot(desired_angle,ILD)
hold on


amp1 = 10^(60/20)*s1_level;
amp2 = 10^(60/20)*s2_level;
amp = amp1.^2+amp2.^2;
dB = 10*log10(amp);

subplot(236)
plot(desired_angle,dB)
hold on
% ylim([-5 5])



%%

