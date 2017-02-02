function pa_plotspec(snd,Fs)
% PA_PLOTSPEC(SND,FS)
% close all;
T = (1:length(snd))/Fs;
 subplot(221);
plot(T,snd,'k-')
ylabel('Amplitude (au)');
xlabel('Time (s)');
xlim([min(T) max(T)]);
hold on
axis square;

subplot(224)
pa_getpower(snd,Fs,'orientation','x');
set(gca,'XTick',[0.5 1 2 4 8 16]*1000);
set(gca,'XTickLabel',[0.5 1 2 4 8 16]);
ax = axis;
xlim(0.6*ax([1 2]));
xlim([100 20000])
set(gca,'XScale','log');
axis square;

subplot(223);
nsamples	= length(snd);
t			= nsamples/Fs*1000;
dt			= 5;
nseg		= t/dt;
segsamples	= round(nsamples/nseg); % 12.5 ms * 50 kHz = 625 samples
noverlap	= round(0.6*segsamples); % 1/3 overlap
window		= segsamples+noverlap; % window size
nfft		= 1000;
spectrogram(snd,window,noverlap,nfft,Fs,'yaxis');
axis square;
% colorbar
cax = caxis;
caxis([0.7*cax(1) 1.1*cax(2)])
set(gca,'YTick',[0.5 1 2 4 8 16]*1000);
set(gca,'YTickLabel',[0.5 1 2 4 8 16]);
set(gca,'YScale','log');
xlim([min(T) max(T)]);
ylim([100 20000]);
drawnow
xlabel('Time (s)');

% linkaxes(ax,'x');