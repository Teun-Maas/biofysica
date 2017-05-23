function tmp
close all
clearvars

% Generates predicted hemodynamic response
fs		= 10;
Nstim	= 1;
on		= 1:Nstim;
on		= on*20; % sec
off		= on+3; % sec
on = on*fs;
off = off*fs;
dur			= 60; % sec
Nsamples	= dur*fs;
X       = zeros(1, Nsamples);
for ii	= 1:length(on)
	X(on(ii):off(ii)) = 1;
end
hemo = canhemofunction(1,X);
pupil = canpupilfunction(1,X);

hemo = hemo/max(hemo);
pupil = pupil/max(pupil);

t = 1:Nsamples;
t = t/fs;
t = t-20;
whos X t
plot(t,X,'k-')
hold on
plot(t,hemo,'r-');
plot(t,pupil,'b-');
ylim([-0.2 1.2]);
function [Y,X] = canhemofunction(beta,X,varargin)
% Y = PA_NIRS_HDRFUNCTION(BETA,X)
%
%

if nargin<1
	beta = [1 1 1];
end
if nargin<2
	X = 1:1000;
end

Fs = keyval('Fs',varargin,10);

%% Initilization
% Fs      = 10;
% gain    = beta(1);
% shift   = 0; % 2 sec response delay
% A		= 4.30; % peak time 4.5 s after stim
% B		= .75;
% N       = length(X); % samples
% dur     = 50; % sec

%%
% x       = 0:(1/Fs):dur; % sec
% x       = x - shift;
% y       = gampdf(x,A,B);
% Y       = conv(X,y);
% Y       = Y(1:N);
% Y       = gain.*Y/max(Y); % Amplitude


%% Hemodynamic impulse response
gain    = beta(1);
if numel(beta)>1
b1 = beta(2);
b2 = beta(3);
else
b1 = 1;
b2 = 1;
end
N       = length(X); % samples

dur     = 50; % sec
xh       = 0:(1/Fs):dur; % sec
a1		= 6; % peak time 4.5 s after stim, shape
a2		= 16; % peak time 4.5 s after stim, shape
c		= 1/6;
y       = gampdf(xh,a1,b1)-c*gampdf(xh,a2,b2);
Y       = conv(X,y);
Y       = Y(1:N);
Y       = gain.*Y; % Amplitude

function [Y,X] = canpupilfunction(beta,X,varargin)



if nargin<1
	beta = [1 1 1];
end
if nargin<2
	X = 1:1000;
end

Fs = keyval('Fs',varargin,10);



%% Hemodynamic impulse response
gain    = beta(1);
if numel(beta)>1
b1 = beta(2);
b2 = beta(3);
else
b1 = 1;
b2 = 1;
end
N       = length(X); % samples

dur     = 50; % sec
xh       = 0:(1/Fs):dur; % sec
a1		= 6; % peak time 4.5 s after stim, shape
a2		= 16; % peak time 4.5 s after stim, shape
c		= 1/6;
y       = gampdf(xh,a1,b1)-c*gampdf(xh,a2,b2);
whos y
w = 10.1;
tmax = 930/1000;
y = 1e4*(xh.^w).*exp(-xh.*w/tmax);
whos X y

Y       = conv(X,y);
Y       = Y(1:N);
Y       = gain.*Y; % Amplitude
