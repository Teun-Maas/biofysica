function movingsound_rationale(spk1,spk2,panningmeth)

close all;
clearvars;

%% Speaker locations
if nargin<1
	spk1			= 0;
end
if nargin<2
	spk2			= 5;
end
if nargin<3
	panningmeth		= 2;
end

%%
desired_angle	= spk1:0.1:spk2; % in deg, at a 'speed' of 0.1 deg/sample
desired_angle	= desired_angle(1:end-1);
s1_lin			= (desired_angle-min(desired_angle))/(max(desired_angle)-min(desired_angle));

%% Head shadow
HSfun		= @(x) 9.7*sin(0.02*x+0.27);
HS			= HSfun(desired_angle);

%% establish loudspeaker locations and numbers
% determine panning levels that ensure a constant overall level
% scaling factors are for scaling amplitude (e.g. in rmsPa)
if panningmeth == 1
	% method 1, used by Brimojoin: sin2+cos2 = 1
	dx				= abs(spk2-spk1); % neighbour speaker separation
	s1_level	= sin(abs(rem(desired_angle,dx)/dx*(pi/2)));
	s2_level	= cos(abs(rem(desired_angle,dx)/dx*(pi/2)));
elseif panningmeth == 2
	% Method 2, thought up by Lof, sqrt(x)+sqrt(1-x) = 1
	s2_lin		= 1-s1_lin;
	s1_level	= sqrt(s1_lin);
	s2_level	= sqrt(s2_lin);
end

spk1lvl		= 60; % in dB
spk2lvl		= 60; % in dB

HSspk1		= HS(1)+spk1lvl;
HSspk2		= HS(end)+spk2lvl;


%% ILD
% Acoustic power of two incoherent noises is summed
pow1		= (10^(HSspk1/20)*s1_level).^2+(10^(HSspk2/20)*s2_level).^2;
pow2		= (10^(HSspk1/20)*s2_level).^2+(10^(HSspk2/20)*s1_level).^2;
dB1			= 10*log10(pow1);
dB2			= 10*log10(pow2);
ILD			= dB2-dB1;

%% "Overall" level
pow1		= 10^(spk1lvl/20)*s1_level; % is this correct?
pow2		= 10^(spk2lvl/20)*s2_level;
pow			= pow1.^2+pow2.^2;
dB			= 10*log10(pow);

%% Graphics
figure(1)
clf
subplot(221)
plot(desired_angle,HS,'o-','MarkerFaceColor','w');
xlabel('Location (deg)');
ylabel('Head shadow (dB)');

subplot(222)
plot(s1_lin+1,s1_level,'o-','MarkerFaceColor','w'); % plot as a function of s1_lin: so as a function of speaker number
hold on
plot(s1_lin+1,s2_level,'o-','MarkerFaceColor','w');
xlim([0.8 2.2]);
set(gca,'XTick',[1 2]);
ylabel('amplitude scaling factor');
xlabel('Speaker number');


subplot(223)
plot(s1_lin+1,ILD,'o-','MarkerFaceColor','w')
hold on
xlim([0.8 2.2]);
ylim([-1 1]);
set(gca,'XTick',[1 2]);
ylabel('ILD (dB)');
xlabel('Speaker number');

subplot(224)
plot(s1_lin+1,dB,'LineWidth',2)
hold on
xlim([0.8 2.2]);
% ylim([-1 1]);
set(gca,'XTick',[1 2]);
ylabel('ILD (dB)');
xlabel('Level (dB)');

if exist('nicegraph','file')
	for ii = 1:4
		subplot(2,2,ii)
		nicegraph;
	end
else
	warning('Biofysica Toolbox not in use!');
end
	
