function [predsum]=ic_prediction2(strf,stim,dt,df,f0)

% function [predsum]=prediction2(strf,stim,dt,df,f0);
% strf: strf of cell
% stim: spectrogram of stim (size(stim,1)=size(strf,1);
% dt: time bin	(ms) default: 12.5 ms
% df: frequency bin (octave) (default: 0.25)
% f0: lowest frequency (Hz) default: 250
%

if (nargin < 5)
   f0=250;
end;

if (nargin < 4)
   df=0.25;
end;

if (nargin < 3)
   dt=12.5;
end;

% prediction: pred=SUMf(strf-CONVt-stim);
xpts=size(strf,1); %number of frequency points;
for i=1:xpts
   M(i,:)=conv(strf(i,:),stim(i,:));
end
pred=M;
predsum=mean(M,1);

%plots
Ntstim=size(stim,2);
Ntpred=size(M,2);
Ntstrf=size(strf,2);
Nf=size(strf,1);
xstrf=Ntstrf/Ntpred;
xstim=Ntstim/Ntpred;

tpred=0:dt:dt*(Ntpred-1);
tstrf=0:dt:dt*(Ntstrf-1);
tstim=0:dt:dt*(Ntstim-1);
f1=[0:Nf-1]*df;
fticks=[0:(Nf-1)/4]*df*4;
flabels=round(f0*2.^fticks);

mxstrf=max(max(abs(strf)));
mxstim=max(max(abs(stim)));
mxpred=max(max(abs(pred)));
mxpred2=max(abs(predsum));


%if nargout==0
figure;
clf;

subplot(4,1,1);
imagesc(tstrf,f1,strf,[-mxstrf,mxstrf]);
title(' strf');
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',fticks,'YTickLabel',flabels);
pos=get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*xstrf pos(4)]);
ylabel('frequency (Hz)');

subplot(4,1,2);
imagesc(tstim,f1,stim,[-mxstim,mxstim]);
title('stim');
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',fticks,'YTickLabel',flabels);
pos=get(gca,'position');
set(gca,'position',[pos(1:2) pos(3)*xstim pos(4)]);
ylabel('frequency (Hz)');
colorbar;

subplot(4,1,3);
imagesc(tpred,f1,M,[-mxpred,mxpred]);
title('prediction');
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',fticks,'YTickLabel',flabels);
ylabel('frequency (Hz)');

subplot(4,1,4);
plot(tpred,predsum,'k');
set(gca,'YLim',[-mxpred2,mxpred2],'Xgrid','On');
xlabel('time (ms)');

set(gcf,'Position',[740,90,400,700]);
