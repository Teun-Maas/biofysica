function [resp, SponR] = ic_callscan(DatFile, tbin, Twindow, MaxRatio, Sigmarr)

% FUNCTION resp = callscan(DatFile,<tbin>,<Twindow>, <MaxRatio>)
% March 2001
%  hv


% Huib Versnel/John van Opstal/Marcel Zwiers
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

if nargin<5
   Sigmarr=[0 2.^(0:5)];
end;

if (nargin < 4)
  MaxRatio = 2.0;
end;

if (nargin < 3)
  Twindow = [251:1900];
end;

if (nargin < 2)
  tbin = 5;
end;

[Stim Spk] = ic_apestim(DatFile);

lspk=size(Spk,2);
Nf = size(Spk,1);
nSigma=length(Sigmarr);
hSigma=round(nSigma/2);

% Zoek het auditieve target
for i = 1:Stim(1,2)
  if Stim(1,i*9-6)==2, snr = i*9+2; end
end

[F,Ind]    = sort(Stim(:,snr));

PreResp   = 1000 * mean(Spk(Ind, 1:Twindow(1)-1)');
PostResp  = 1000 * mean(Spk(Ind, Twindow(end)+1:lspk)');
Resp      = 1000 * mean(Spk(Ind, Twindow)');

Spks=Spk(Ind,:);
%PreResp
SponR = mean(PreResp);

MinResp  = 50;   % spikes/s

for i = 1:Nf

  inear = neighbor(i, 1, Nf);
  PreRatio(i)  = (PreResp(i)  + MinResp)/(MinResp + SponR);
  Ratio(i)     = (Resp(i)     + MinResp)/(MinResp + mean(Resp(inear)));
  PostRatio(i) = (PostResp(i) + MinResp)/(MinResp + mean(PostResp(inear)));

  IntraOK = Ratio(i)    < MaxRatio;
  PreOK   = PreRatio(i) < MaxRatio;
  PostOK  = PreRatio(i) < MaxRatio;

  if (~IntraOK || ~PreOK || ~PostOK)
     Resp(i) = min(Resp(inear(inear~=i)));
     disp (['  Correction at F = ', int2str(F(i)) ' Hz']);
  end;

end;


jspl=mod(F-90000,10); spl=10*(7-jspl);
Fmonk=F(F<90100);  Spkmonk=Spks(F<90100,:); splmonk=spl(F<90100);
Fbird=F(F>90100);  Spkbird=Spks(F>90100,:); splbird=spl(F>90100);
splmonk=splmonk(1:5:length(splmonk));
magfac=1000;

for k=1:6
   q=nanmean(Spks(k*15-14:k*15,:),1)*magfac;
   psths(k,:)=ic_rebinfrac(q,tbin);
end
ypsths=10*ceil(max(max(psths))/10);

for k=1:18,
   q=mean(Spks(k*5-4:k*5,:),1)*magfac;
   psthlevs(k,:)=ic_rebinfrac(q(Twindow),tbin);
end

ypsthlevs=10*ceil(max(max(psthlevs))/10);

for k=1:3,
   labely(k,:)=[num2str(splmonk(k)) ' dB'];
end

%plots
nticks=floor(length(Fmonk)/3);
timetick=[0:250:1750];

% titles(1,:)=[DatFile, '      Grunt'];
% titles(2,:)=['     Pulsed'];
% titles(3,:)=[' Undulating'];
% titles(4,:)=[DatFile, '  Blackbird'];
% titles(5,:)=[' Meadowlark'];
% titles(6,:)=['     Oriole'];

Twind{1}=[275:562]; Twind{2}=[275:600]; Twind{3}=[275:1534];
Twind{4}=[275:1484]; Twind{5}=[275:1750]; Twind{6}=[275:1894];

if nargout==0
figure(1); clf
for k=1:3,
   subplot(3,1,k);
   ispk=(1:15)+15*(k-1);
   dotplot(Spks(ispk,:),0);
%    title(titles(k,:));
   set(gca,'ytick', [1:5:nticks], 'TickDir','out','FontSize',10);
   set(gca,'yticklabel', labely,'XGrid','off');
   set(gca,'XTick', timetick,'XTicklabel',num2str(timetick'/1000));
end;

Xlabel('Time (ms)');
set(gcf,'position',[300 10 700 750]);

for v=1:6,
	for m=1:nSigma,
   	disp(['Sigma: ',num2str(Sigmarr(m))]);
   	Spd=ic_spkdenshv(Spks((1:15)+(v-1)*15,Twind{v}),Sigmarr(m));
   	cc=corrcoef(Spd');
      for k=1:3,
         j=(1:5)+5*(k-1);
%         cau(k+(v-1)*3,m)=meancorrelation(cc(j,j));
          cau(m,v,k)=ic_meancorrelation(cc(j,j));
      end;
   end;
end;
%cau(hSigma,:,:)
for k=1:3,
   [mx,p]=max(mean(cau(:,:,k),2));
   maxSigma=Sigmarr(p);
   maxC(k)=mx;
   maxSig(k)=maxSigma;
end;
disp('maximum correlations averaged over all calls');
maxC';
disp('best Sigma (ms)')
maxSig';
Sigmah=Sigmarr(hSigma);
disp(['correlations at highest level, at ',num2str(Sigmah),' ms']);
cau(hSigma,:,1)';

figure(2); clf
for k=1:3,
   subplot(3,1,k);
   ispk=(1:15)+15*(k-1);
   dotplot(Spkbird(ispk,:),0);
   title(titles(k+3,:));
   set(gca,'ytick', [1:5:nticks], 'TickDir','out','FontSize',10);
   set(gca,'yticklabel', labely,'XGrid','off');
   set(gca,'XTick', timetick,'XTicklabel',num2str(timetick'/1000),'TickDir','out');
end;

Xlabel('Time (ms)');
set(gcf,'position',[300 10 700 750]);

figure(3); clf

lenp=size(psths,2);   % length of psth
timetick=[0:lenp/4:lenp];

for k=1:6,
   m=mod(k+2,3).*2+ceil(k/3);
   subplot(3,2,m);
   bar(psths(k,:),'k');
   title(titles(k,:));
   ax=axis;
	axis([ax(1) lenp ax(3) ypsths]);
   set(gca,'XTick', timetick,'XTicklabel',num2str(tbin*timetick'/1000),'TickDir','out');
end

set(gcf,'position',[300 10 700 750]);
  
figure(4); clf

lenp=size(psthlevs,2);   % length of psth
timetick=[0:lenp/4:lenp];

for k=1:18,
   m=mod(k+8,9).*2+ceil(k/9);
   subplot(9,2,m);
   bar(psthlevs(k,:),'k');
   ax=axis;
	axis([ax(1) lenp ax(3) ypsthlevs]);
   m=div(k+2,3);
	if mod(k,3)==1 
      title(titles(m,:));
   end;
	if mod(k,9)==0
      set(gca,'XTick', timetick,'XTicklabel',num2str(tbin*timetick'/1000),'TickDir','out');
      Xlabel('Time (ms)');
   else
      set(gca,'XTick', timetick,'XTicklabel',' ','TickDir','out');
   end  
   m=div(mod(k-1,9),3);
   hpos=get(gca,'position');
	set(gca,'position',hpos+[0,-.01*m,0,0]);
end

set(gcf,'position',[300 10 700 750]);

figure(5)
clf;
for v=1:6,
   m=2*v-1-5*div(v,4);
   subplot(3,2,m);
	%m=3*(v-1);
	plot(Sigmarr,cau(:,v,1),'g-',Sigmarr,cau(:,v,2),'m-',Sigmarr,cau(:,v,3),'c-','linewidth',1.5);
	xlabel('Sigma (ms)');
   ylabel('intertrial correlation');
end;
set(gcf,'position',[300 10 700 750]);


else 
   for k=1:6,
      for l=1:3,
      	m=(k-1)*3+l;
      	resp{k,l}=psthlevs(m,:);
      end;
   end;
end;


     

% -------- plot a line indicating background activity (prior to stimulus);

%BGlevel   = 1000 * mean(mean(Spk(:,1:200)));
%Meanlevel = 1000 * mean(mean(Spk(:,Twindow)));

% -------- plot the bar graph

%plot([1:length(Resp)], Resp);
%box off

%horzlijn(BGlevel, 'r:');

%set(gca,'xtick', [1:3:length(Resp)], 'TickDir','out');
%set(gca,'xticklabel', F(1:3:end));

%xlabel('Frequency (Hz)');
%ylabel('Mean firing rate (Hz)');

%set(gcf, 'numbertitle', 'off');
%set(gcf, 'name', DatFile);

%[M,I] = max(Resp);
%Fmax  = F(I);
%disp([sprintf('Fmax at: %d Hz',Fmax)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ivec = neighbor (i, istart, istop)
%
% Return the neighborhood of i.
% For i that is [i-1 i i+1] normally, only at the edge
% when eg i=1 it will be [1 2 3].

if (i<=istart)
  ivec = istart + [0 1 2];
elseif (i>=istop)
  ivec = istop - [2 1 0];
else
  ivec = [i-1 i i+1];
end;

