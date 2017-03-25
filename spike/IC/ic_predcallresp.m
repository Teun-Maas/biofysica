function ic_predcallresp(strf,DatFile,spcs,dt,df,f0,loglin)
%
%function predcallresp(strf,DatFile,spcs,dt,df,f0,loglin);
% strf: from strf1 or strf2
% Datfile: file with responses to calls
% spcs: spectra of calls, frequencies match those of strf
% 		structure s{k} with k=1:6, obtained with specvocalsalt
% dt: time bin (default: 12.5 ms)
% df : frequency bin (default: 0.25 oct);
% f0: lowest frequency (plotting purpose)
% loglin: 0 for lin, 1 for log (default);
%

% Huib Versnel/John van Opstal/Marcel Zwiers
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<7
   loglin=1;
end;
if nargin<6
   f0=250;		%plotting purpose
end;
if nargin<5
   df=0.25;
end;
if nargin<4
   dt=12.5;
end;
[pth,nm]	= fileparts(DatFile);
disp(DatFile);
flpar		= fullfile(pth,['paramcall_',nm ,'.txt']);
flpredresp	= fullfile(pth,['predrespcall_',nm ,'.txt']);

%% Spike data
[resp, sr] = ic_callscan(DatFile, dt);  %resp{k} with k=1:6 representing 6 stimuli

titles(1,:) = [DatFile, '      Grunt'];
titles(2,:) = [DatFile, '     Pulsed'];
titles(3,:) = [DatFile, ' Undulating'];
titles(4,:) = [DatFile, '  Blackbird'];
titles(5,:) = [DatFile, ' Meadowlark'];
titles(6,:) = [DatFile, '     Oriole'];

prd			= cell(6,1);
predrate	= NaN(6,1);
for k = 1:6,
   spec1 = spcs{k};
   if loglin,
      spec1 = log10(spec1+1);	% add 1 to prevent negative numbers
   end;
   prd{k}			= ic_prediction2(strf,spec1,dt,df,f0);
   predrate(k)		= 1000*mean(prd{k}(prd{k}>0))/dt;
   for j = 1:3,
       res				= resp{k,j};					
       n				= min(length(prd{k}),length(res));
       resp1			= res(1:n)-sr;
       resprate1(k,j)	= mean(resp1(resp1>0));
       resprate2(k,j)	= mean(res(1:n));
       rmsresp(k,j)		= sqrt(mean(resp1.^2));
   end
   rmspred(k)=sqrt(mean(prd{k}.^2));
%    pause;
end
for k=1:6,
   figure;
   disp(['rms response: ',num2str(round(rmsresp(k)))]);
   p=prd{k};
   p1=p; p1(p1<0)=0;
   C=plotpredcall(resp,p,k,titles,sr,dt);
   Cs(:,k)=C';
   figure;
   C=plotpredcall(resp,p1,k,titles,0,dt);
   Csrect(:,k)=C';
%    pause;
end;

for j=1:3,
    c1=corrcoef(predrate,resprate1(:,j));
    c2=corrcoef(predrate,resprate2(:,j));
    c3=corrcoef(rmsresp(:,j),rmspred);
    r1(j)=c1(1,2);
    r2(j)=c2(1,2);
    r3(j)=c3(1,2); 
end

[~,kr1]=max(resprate1(:,1));
[~,kr3]=max(rmsresp(:,1));
[~,kp1]=max(predrate(:,1));
[~,kp3]=max(rmspred);

disp('rms response: ');
for k=1:6,
   disp(num2str(rmsresp(k,1)));
end
disp(' ');
disp([' max of mean response>bg: call',num2str(kr1)]);
disp([' max of rms response: call',num2str(kr3)]);
disp([' max of prediction>0: call',num2str(kp1)]);
disp([' max of rms prediction: call',num2str(kp3)]);
disp(' ');

p=polyfit(resprate1(:,1),predrate,1);
y1=p(1)*resprate1(:,1)+p(2);

figure;
subplot(2,1,1);
plot(resprate1,predrate,'*',resprate1(:,1),y1,'-');
xlabel('mean response > background');
ylabel('mean prediction > 0');
xcord=get(gca,'XLim'); x=(xcord(2)-xcord(1))*0.2+xcord(1);
ycord=get(gca,'YLim'); y=(ycord(2)-ycord(1))*0.8+ycord(1);
text(x,y,['r: ',num2str(r1)]);
disp(['mean response>bg vs mean prediction>0:  ',num2str(r1)]);
disp(['mean response vs mean prediction>0:  ',num2str(r2)]);

p=polyfit(rmsresp(:,1)',rmspred,1);
y2=p(1)*rmsresp(:,1)+p(2);
   
subplot(2,1,2);
plot(rmsresp,rmspred,'*',rmsresp(:,1),y2,'-');
xlabel('rms response-background');
ylabel('rms prediction');
xcord=get(gca,'XLim'); x=(xcord(2)-xcord(1))*0.2+xcord(1);
ycord=get(gca,'YLim'); y=(ycord(2)-ycord(1))*0.8+ycord(1);
text(x,y,['r: ',num2str(r3(1))]);
disp(['rms response-bg vs rms prediction:  ',num2str(r3(1))]);

par=[Cs Csrect rmsresp' [rmspred; rmspred; rmspred] [r1; r2; r3]];
predresp=[rmsresp rmspred'];

save(flpar,'par', '-ascii');
save(flpredresp,'predresp', '-ascii');

disp(' ');
disp(['correlation coefficients saved to ',flpar]);
disp(['rms response and prediction saved to ',flpredresp]);


function C=plotpredcall(histf,pred,kcall,titles,sr,dt)

% function C=plotpredcall(histf,pred,kcall,sr,dth,dtp);
%
%   C: correlation between prediction and response
%

if (nargin < 6)
   dt=12.5;
end;

if (nargin < 5)
   sr=0;
end;

if (nargin < 4)
   titles=repmat([' '],6);
end;

for l=1:3,
   r(l,:)=histf{kcall,l}-sr;
end;

nr=size(r,2);
np=length(pred);
yscl=max(max(r));

pred1=yscl/max(pred)*pred;
n=min(np,nr);

for l=1:3,
   r1(l,:)=r(l,1:n);
end;
pred1=pred1(1:n);
tr=0:dt:(n-1)*dt;

scl1=yscl*1.1;


for l=1:3,
    cc=corrcoef(r1(l,:),pred1);
    C(l)=cc(1,2);   
end;


disp(['correlations for call',num2str(kcall),': ']);
for k=1:3,
   disp(num2str(C(k)));
end;
disp(' ');

subplot(3,1,1);
plot(tr,r1(1,:),'m',tr,pred1,'c',tr,0*tr,'k--','linewidth',1.5);
set(gca,'YLim',[-scl1,scl1],'TickDir','out','Xticklabel',[],'Box','off');
title(titles(kcall,:));

subplot(3,1,2);
plot(tr,r1(2,:),'m',tr,pred1,'c',tr,0*tr,'k--','linewidth',1.5);
set(gca,'YLim',[-scl1,scl1],'TickDir','out','Xticklabel',[],'Box','off');

subplot(3,1,3);
plot(tr,r1(3,:),'m',tr,pred1,'c',tr,0*tr,'k--','linewidth',1.5);
set(gca,'YLim',[-scl1,scl1],'TickDir','out','Box','off');
set(gcf,'position',[500,35,650,760]);
xlabel('time (ms)');
