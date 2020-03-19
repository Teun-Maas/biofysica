function phon = spl2phon(spl,f)
%
% PHON = SPL2PHON(SPL,F)
%
% Determines phon level from a pure tone of frequency F presented at a
% level of SPL dB(SPL).
%
% Generates phon levels from Equal Loudness Contours, interpolated from
% values as described in ISO226. 
%
% The input frequencies F can range from 20Hz - 12.5kHz.
%
% The valid SPL range is from 0 - 90 dB SPL. Values outside this range do
% not have experimental values. 
%
% Based on iso226.m from Jeff Tackett 03/01/05
% https://nl.mathworks.com/matlabcentral/fileexchange/7028-iso-226-equal-loudness-level-contour-signal
%
% See also ISO226, GRIDDEDINTERPOLANT

%%
if nargin<1
	spl = 60;
end
if nargin<2
	f = 5000;
end

fname = 'spl2phon.mat';
if exist(fname,'file')
	load(fname,'F');
else
	%%  ISO226 Tables
	freq = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 ...
		1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500];
	
	af = [0.532 0.506 0.480 0.455 0.432 0.409 0.387 0.367 0.349 0.330 0.315 ...
		0.301 0.288 0.276 0.267 0.259 0.253 0.250 0.246 0.244 0.243 0.243 ...
		0.243 0.242 0.242 0.245 0.254 0.271 0.301];
	
	Lu = [-31.6 -27.2 -23.0 -19.1 -15.9 -13.0 -10.3 -8.1 -6.2 -4.5 -3.1 ...
		-2.0  -1.1  -0.4   0.0   0.3   0.5   0.0 -2.7 -4.1 -1.0  1.7 ...
		2.5   1.2  -2.1  -7.1 -11.2 -10.7  -3.1];
	
	Tf = [ 78.5  68.7  59.5  51.1  44.0  37.5  31.5  26.5  22.1  17.9  14.4 ...
		11.4   8.6   6.2   4.4   3.0   2.2   2.4   3.5   1.7  -1.3  -4.2 ...
		-6.0  -5.4  -1.5   6.0  12.6  13.9  12.3];
	
	%% Create Grid
	x			= 1:90;
	[~,af]		= ndgrid(x,af);
	[~,Lu]		= ndgrid(x,Lu);
	[~,Tf]		= ndgrid(x,Tf);
	[x,y]	= ndgrid(x,freq);
	
	%% SPL2Phon
	Af			= 10.^(af.*(x+Lu-94)/10);
	v			= 40*log10( (Af-(0.4*10.^(((Tf+Lu)/10)-9 )).^af)/4.47e-3 + 1.15 );
	
	%% Grid Interpolation
	F			= griddedInterpolant(x,y,v,'spline');
	vq			= F(spl,f);
	
% 	save('/Users/marcw/Gitlab/biofysica/audition/spl2phon.mat','F');
	save('spl2phon.mat','F');
end

%% Output
phon = F(spl,f);

%% Some graph
% close all
% xq		= 1:90;
% yq		= 100:10000;
% [xq,yq] = ndgrid(xq,yq);
% vq		= F(xq,yq);
% plot(yq',vq');

