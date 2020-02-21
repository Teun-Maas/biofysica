%% Cochlear filter bank
%
% We will build two cochleas, each with 1000 hair cells, with a tuning curve
% described by a gamma-tone filter each with their own characteristic
% frequency F_C.
%
% ILDs will be processed by taking the difference in response
%
% Cross-channel integration will take place before or after this
close all;
clearvars;

figure(1)

Fc			= oct2bw(1000,0:0.1:4); % characteristic frequency
nUnits		= numel(Fc);
Fs			= 40000;
nsamples	= 2^12;
fcoefs		= MakeERBFilters(40000,Fc);
y			= ERBFilterBank([1 zeros(1,nsamples-1)], fcoefs);
resp		= abs(fft(y'));
freq		= (0:nsamples-1)/nsamples*Fs;
freq		= freq(1:nsamples/2-1)';
resp		= resp(1:nsamples/2-1,:);
semilogx(freq,resp);
% axis([100 16000 -60 0])
xlabel('Frequency (kHz)');
ylabel('Filter Response (dB)');
nicegraph;
axis normal
set(gca,'XTick',[1000 2000 4000 8000 16000],'XTickLabel',[1000 2000 4000 8000 16000]/1000);
xlim([500 25000]);




%% What else
%
% So?
%

figure(2)

freq	= (0:511)/512*40000;
freq	= freq(1:255)';
nfreq	= numel(freq);
Z		= eye(nfreq,nfreq);

Level			= 0:100;
Pres		= 20*10.^(Level/20);
nPres		= numel(Pres);

nUnit = 16;
fcoefs	= MakeERBFilters(40000,nUnit,2000);
P = NaN(nUnit,nPres,nfreq);
s = 1;
th = 0;
for jj = 1:nUnit
	for ii = 1:nPres
		y		= ERBFilterBank([Pres(ii) zeros(1,511)], fcoefs);
		resp	= abs(fft(y'));
		resp	= resp(1:255,:);
		f		= max(Z.*resp(:,jj));
		p		= s*log(1+exp(log10(f)-th)/s);
		
		%%
% 		hold on
% 		plot(f,p)
		
% 		q		= p-0.25*p;
		P(jj,ii,:) = p;
	end
		%%
% 		keyboard
	figure(2)
	subplot(4,4,jj)
	hold on
	nicegraph;
	cmap = (cbrewer('seq','Reds',64,'pchip'));
	colormap(cmap);
	contourf(freq,Level,squeeze(P(jj,:,:)),'k');
	set(gca,'Xscale','log','XTick',[1000 2000 4000 8000 16000],'XTickLabel',[1000 2000 4000 8000 16000]/1000);
	[mx,idx] = max(max(squeeze(P(jj,:,:))));
	Fc(jj) = freq(idx);
	xlim([500 20000]);
	
% 	figure(3)
% 	hold on
% 	whos freq P
% 		plot(Level,squeeze(P(jj,:,125)),'k');

end

%% ILD dependence
% Note cells in LSO are predominantly EI
ild = -10:2:10;
nild = numel(ild);
	figure(5)
clf
hold on
base = 60;
for ii = 1:nild
	sel		= Level==base+ild(ii);
	Tipsi	= squeeze(P(:,sel,:));
	sel		= Level==base-ild(ii);
	Tcontra = squeeze(P(:,sel,:));
	T = Tipsi-0.25*Tcontra;
	
	for jj = 1:nUnit
% 	figure(5)
% 	subplot(4,4,jj)
% 	hold on
% 	semilogx(freq,T(jj,:))
% 	set(gca,'Xscale','log','XTick',[1000 2000 4000 8000 16000],'XTickLabel',[1000 2000 4000 8000 16000]/1000);
% 	xlim([500 20000]);
	[mx,idx] = max(T(jj,:));
	ILD(jj,ii) = mx; 
	end
	

end

%%
figure(6)
clf
plot(ild,ILD','o-')
nicegraph
%%
return
%%
%% Cochlear filter bank
%
% We will build two cochleas, each with 1000 hair cells, with a tuning curve
% described by a gamma-tone filter each with their own characteristic
% frequency F_C.
%
% ILDs will be processed by taking the difference in response
%
% Cross-channel integration will take place before or after this
clearvars;
Level		= 0:100; % level in contralateral ear
Pres		= 20*10.^(Level/20);
nPres		= numel(Pres);
