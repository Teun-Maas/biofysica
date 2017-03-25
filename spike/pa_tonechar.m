function ToneChar = pa_tonechar(fname,varargin)
% PA_TONECHAR(FNAME,METH,DSP)
%
% Characterize tone response characteristics of auditory neurons.
%
% Characteristics are as defined by Recanzone, Guard and Chan, 2000, and
% include:
%
% - Frequency response are
% - spontaneous rate
% - driven response
% - characteristic frequency


% (c) 2012 Marc van Wanrooij

%% Initialization
dsp = pa_keyval('display',varargin);
if isempty(dsp)
	dsp = 1;
end
sd = pa_keyval('sd',varargin);
if isempty(sd)
	sd = 2;
end
% if nargin<2
% 	meth = 2;
% end

%% Load f32-data set
dat			= spikematf(fname,1);
toneonset	= 300; % Hard-coded. Where can this be found in stimulus file?
charDur		= 50; % time (ms) during stimulus presentation, used for
% characterization. Ideally, this should equal sound duration. However, for
% onset response cells (as opposed to sustained response cells), this might
% skew the characterization.
preDur		= 300;

%% Obtain firing rate
N			= NaN(1,length(dat)); % number of spikes
Frequency	= N; % frequency
soundLevel	= N; % sound level (dB)
soundDur	= N; % sound duration (ms)
B			= N;
firingRate	= N; % Firing Rate during first 50 ms (spikes/s)
stimRept	= N; % number of repeats
A			= [];
spike		= struct([]);
for ii = 1:length(dat)
	Frequency(ii)	= dat(ii).stim(1);
	soundLevel(ii)	= dat(ii).stim(2);
	soundDur(ii)	= dat(ii).stim(3);
	
	b				= full(mean(dat(ii).sweeps));
	a				= full([dat(ii).sweeps]);
	A				= [A;a]; %#ok<AGROW>
	t				= 1:length(b);
	
	sel				= t>=toneonset & t<toneonset+charDur;
	N(ii)			= sum(b(sel)); % number of spikes
	stimRept(ii)	= size(a,1);
	B(ii)			= sum(b(~sel));
	firingRate(ii)	= N(ii)/charDur*1000; % Firing Rate during first 50 ms (spikes/s)
	b				= a(:,sel);
	spike(ii).spiketime	= b';
	spike(ii).Frequency = repmat( dat(ii).stim(1),1,size(b,1));
	spike(ii).instantfiringrate = mean(b)'/charDur*1000;
end
uFrequency		= unique(Frequency);
nFrequency		= numel(uFrequency);
uLevel			= unique(soundLevel);
nLevel			= numel(uLevel);
rept			= unique(stimRept);
if numel(rept)>1
	rept = mean(rept);
	disp('Number of repeats vary');
end

%% Spike density
SDF = pa_spk_sdf(A);
%% Some visual check
if dsp
	figure(101)
	subplot(211)
	pa_spk_dotplot(A,'markersize',5);
	pa_verline([300 450],'r-');
	pa_horline(cumsum(rept*nLevel*ones(nFrequency,1)),'r-');
	
	AF	= [spike.Frequency];
	x	= 1:length(AF);
	set(gca,'YTick',x(1:rept*nLevel:end),'YTickLabel',round(AF(1:rept*nLevel:end)));
	plot(rept*nLevel*nFrequency*SDF/max(SDF),'r-','LineWidth',2);
end

%% Interpolate
% intFlag = true;
% if intFlag
% 	x	= soundLevel;
% 	y	= Frequency;
% 	z	= firingRate;
% 	x	= reshape(x,nLevel,nFrequency)';
% 	y	= reshape(y,nLevel,nFrequency)';
% 	z	= reshape(z,nLevel,nFrequency)';
% 	XI	= linspace(min(soundLevel),max(soundLevel),7);
% 	YI	= log2(pa_oct2bw(250,0:0.25:6));
% 	[XI,YI] = meshgrid(XI,YI);
% 	ZI	= interp2(x,log2(y),z,XI,YI,'linear');
% 
% 
% 	keyboard
% 	x	= soundLevel;
% 	y	= Frequency;
% 	z	= [spike.instantfiringrate]';
% 	x	= reshape(x,nLevel,nFrequency)';
% 	y	= reshape(y,nLevel,nFrequency)';
% 	n = size(z,2);
% 	x = repmat(x,[1,1,n]);
% 	y = repmat(y,[1,1,n]);
% 	t = repmat(1:n,[nFrequency,1,nLevel]);
% 	t = permute(t,[1,3,2]);
% 	z	= reshape(z,[nFrequency,nLevel,n]);
% 	XIi	= linspace(min(soundLevel),max(soundLevel),7);
% 	YIi	= log2(pa_oct2bw(250,0:0.25:6));
% 	TI = 1:n;
% 	[XIi,YIi,TI] = meshgrid(XIi,YIi,TI);
% 	II	= interp3(x,log2(y),t,z,XIi,YIi,TI,'linear');
%    
% 	firingRate = ZI(:);
% 	soundLevel = XI(:);
% 	Frequency = 2.^YI(:);
% 	uFrequency		= unique(Frequency);
% 	nFrequency		= numel(uFrequency);
% 	uLevel			= unique(soundLevel);
% 	nLevel			= numel(uLevel);
% 	instFiringRate = reshape(II,nFrequency*nLevel,n);
% end

z	= firingRate;
z	= reshape(z,nLevel,nFrequency);
sz = medfilt2(z,[2 1]);
% h = fspecial('average', [2 1]);
% sz = imfilter(h,z);
% subplot(122);imagesc(sz)
% subplot(121);imagesc(z)
% firingRate = reshape(sz,1,nLevel*nFrequency);
firingRate = sz(:)';

% keyboard

%% Spontaneous activity
% Measured in 50 ms before stimulus onset
t				= 1:size(A,2); % time (ms)
sel				= t>toneonset-preDur & t<toneonset; % prestimulus period
sponAct			= sum(A(:,sel),2)/preDur*1000; % spontaneuous rate (spikes/s) per trial
muSponAct		= mean(sponAct); % mean over trials
sdSponAct		= std(sponAct); % sd over trials
% keyboard
%% Driven response
% average of all presentations of each stimulus that was greater than mean
% and 2 SDs of the spontaneous activity
minDrivenAct	= muSponAct+sd*sdSponAct; %(spikes/s)
drivenResponse	= firingRate > minDrivenAct; % is there a driven response?
drivenLevel		= soundLevel(drivenResponse); % what sound levels drive a response?
drivenFrequency	= Frequency(drivenResponse); % what sound levels drive a response?

%% Check whether the lowest driven level also induces responses at higher levels
udF = unique(drivenFrequency);
checkLF = NaN(numel(udF),2);
for ii = 1:numel(udF)
	sel = drivenFrequency==udF(ii);
	mn = min(drivenLevel(sel));
	if	mn<70
		if sum(sel)==1
			checkLF(ii,1) = mn;
			checkLF(ii,2) = udF(ii);
		end
	end
end
sel = isnan(checkLF(:,1));
checkLF = checkLF(~sel,:);

if ~isempty(checkLF)
	sel1 = ismember(drivenLevel,checkLF(:,1));
	sel2 = ismember(drivenFrequency,checkLF(:,2));
	sel = sel1&sel2;
	drivenLevel = drivenLevel(~sel);
end
% ismember
% checkLF

if dsp
	figure(101)
	subplot(211)
	str = {fname;['Spon A = ' num2str(round(muSponAct)) ' spikes/s'];['Driven A = ' num2str(round(minDrivenAct)) ' spikes/s']};
	title(str)
end

%% Threshold
threshold		= min(drivenLevel); % lowest intensity that drives a response

%% Characteristic Frequency
% the frequency that produced a driven response at the lowest intensity
indx			= find(drivenLevel == threshold);
charFrequency	= Frequency(drivenResponse);
charFiringRate	= firingRate(drivenResponse);
charFrequency	= charFrequency(indx);
charFiringRate	= charFiringRate(indx);
sel				= ismember(uFrequency,charFrequency);
if sum(diff(sel)<0)>1 % multiple non-adjacent driven frequencies
	[~,indx]		= max(charFiringRate);
	charFrequency	= charFrequency(indx);
else
	charFrequency	= sum(charFrequency.*charFiringRate)/sum(charFiringRate); % linearly averaged characteristic firing rate
end

%% Latency
[~,indx]				= min(abs(uFrequency-charFrequency));
closeCharFrequency		= uFrequency(indx); % the frequency actually tested closest to CF
indx					= indx-2:indx+2; % choose 5 frequencies nearest to CF
indx(indx<1)			= 1; % if extreme frequency, limit the 5 nearest frequencies
indx(indx>nFrequency)	= nFrequency;
nCharFrequency			= numel(indx);
latFrequency			= uFrequency(indx);

n				= numel(latFrequency);
A				= [];
F				= [];
for ii = 1:n
	sel = Frequency==latFrequency(ii);
	tmp = [spike(sel).spiketime];
	A	= [A tmp]; %#ok<AGROW>
	tmp = [spike(sel).Frequency];
	F	= [F tmp]; %#ok<AGROW>
end
charSDF = pa_spk_sdf(A');

if dsp
	figure(101)
	subplot(223)
	pa_spk_dotplot((A'),'markersize',5);
	pa_verline(50,'r-');
	pa_horline(cumsum(rept*nLevel*ones(nCharFrequency,1)),'r-');
	plot(rept*nLevel*nCharFrequency*charSDF/max(charSDF),'r-','LineWidth',2);
	
	x = 1:length(F);
	set(gca,'YTick',x(1:rept*nLevel:end),'YTickLabel',round(F(1:rept*nLevel:end)));
	
	stairs(rept*nLevel*nCharFrequency*sum(A,2)./max(sum(A,2))/2,'r-');
end

A	= mean(A,2)*1000; % number of spikes per millisecond
A	= reshape(A,2,length(A)/2);
A	= sum(A)/2;
minDrivenAct	= muSponAct+2*sdSponAct; %(spikes/s)

b	= A>minDrivenAct;

c				= b+[0 b(1:end-1)]+[0 0 b(1:end-2)];
onsetLatency	= 2*(find(c==3,1)-2); % onset latency ms
if isempty(onsetLatency)
	onsetLatency = NaN;
end

if isempty(threshold)
	threshold = NaN;
end
ToneChar.onsetLatency	= onsetLatency;
ToneChar.charFrequency	= charFrequency;
ToneChar.muSponAct		= muSponAct;
ToneChar.minDrivenAct	= minDrivenAct;
ToneChar.threshold		= threshold;

%% Tuning bandwidth
% keyboard

x		= reshape(Frequency,nLevel,nFrequency);
x		= log2(x);
y		= reshape(soundLevel,nLevel,nFrequency);
z		= reshape(firingRate,nLevel,nFrequency);

XI		= pa_oct2bw(250,linspace(0,6,31));
XI		= log2(XI); % octave
YI		= 10:10:70;
[XI,YI] = meshgrid(XI,YI);
ZI		= interp2(x,y,z,XI,YI);


if isnan(threshold)
	threshold = 30;
end
sel		= YI == threshold+10;
z		= ZI(sel);
x		= 2.^XI(sel); % Hz
indx1	= find(z>minDrivenAct,1,'first');
indx2	= find(z>minDrivenAct,1,'last');
% dz		= [0;diff(z)];
% ddz		= [0;diff(dz)];

if isempty(indx1) || isempty(indx2)
	BW10 = NaN;
else
	BW10	= pa_freq2bw(x(indx1),x(indx2)); % oct
end
ToneChar.BW10 = BW10;

if dsp
	figure(101)
	subplot(223)
	pa_verline(onsetLatency,'r--');
	str = {['CF = ' num2str(round(charFrequency)) ' Hz'];['Onset Latency = ' num2str(round(onsetLatency)) ' (ms)']};
	title(str)
	
	figure(101)
	subplot(211)
	h = pa_verline(300+onsetLatency,'r--');
	set(h,'LineWidth',2);
	
	F = [spike.Frequency];
	indx = find(F==closeCharFrequency);
	pa_horline(indx,'r:');
	
	subplot(224)
	h = semilogx(x,z,'k-'); set(h,'LineWidth',2);
	hold on
% 	h = semilogx(x,ddz-min(ddz),'b-'); set(h,'LineWidth',2);
	
	x = reshape(Frequency,nLevel,nFrequency);
	y = reshape(firingRate,nLevel,nFrequency);
	semilogx(mean(x),mean(y),'k:');
	xlim([min(uFrequency) max(uFrequency)]);
	ylim([0 100]);
	pa_horline(minDrivenAct,'r-');
	f = pa_oct2bw(charFrequency,[-0.5 0 0.5]);
	pa_verline(f,'r:');
	str = {['BW10 = ' num2str(BW10,2)]};
	title(str)
	set(gca,'XTick',uFrequency,'XTickLabel',round(uFrequency));
end
