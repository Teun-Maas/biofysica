 (ff=0.0function pa_createzwuis(hrmnc,ff,fname)
% Function CREATERANDACC
%


%% Initialization
close all;
clear all;
pa_datadir;
Fs					= 50000;          % Hz
Duration            = 20;          % sec
Amp                 = 30;
X                   = 0:1/Fs:Duration;



if nargin <1
    hrmnc = 1;
end
if nargin<2
    ff=2;
end
if length(hrmnc)    == 1
    if hrmnc        == 1
        harmonic    = [2, 3, 7, 11, 19];
%         harmonic    = primes(10)/10;
    elseif hrmnc    == 2
        harmonic    = [5, 7, 9, 16, 23];
    elseif hrmnc    == 3
        harmonic    = primes(32);
    elseif hrmnc    == 4
        harmonic    = 1;
    end
else
    harmonic        = hrmnc;
end

if ff               == 1
    fund            = 0.05;
elseif ff           == 2
    fund            = 0.1;
else
    fund            = ff;
end
freq                = fund.*harmonic;
Y                   = zeros(length(freq),length(X));
InvAmp              = 2*pi*freq;

%% Creating Stimulus
phi = rand(length(harmonic),1); % Random phase
nfreq = length(freq);
for indx               = 1:nfreq
%     Y(i,:)          = Amp.*sin(2*pi*freq(i)*X+phi(i)*2*pi)/InvAmp(i);
    Y(indx,:)          = Amp.*sin(2*pi*freq(indx)*X+phi(indx)*2*pi);
end
Y                   = sum(Y)/nfreq;
Y                   = Y(:);
% Y                   = envelope (Y, NEnvelope);
Y                   = Y';

subplot(211)
plot(X,Y,'k-');
xlabel('Time (s)');
ylabel('Location (deg)');
box off
set(gca,'TickDir','out');

[F,A,PH] = pa_getpower(Y,Fs);

subplot(212)
% sel = F<2;
% semilogx(F(sel),A(sel),'ko-');
% xlabel('Frequency');
% ylabel('Location (deg)');
% box off
% set(gca,'TickDir','out');
% pa_verline(freq);

window      = 2^19;
overlap     = 2^1;
[Pyy,f]     = pwelch(Y,window,overlap,[],Fs);
Pyy         = sqrt(Pyy*Fs);
plot(f,Pyy,'ko-','MarkerFaceColor','w');

ylabel('Power (units)');
xlim([0 3])
pa_verline(freq);

return
Y                   = tableformat(Y,SampleFreq,nrntpulse);

%% Saving
if nargin<3 
    fname     = ['ss' num2str(round(Duration),'%d') 's' num2str(Amp,'%d') 'd' num2str(fund*1000,'%d') 'mhz' num2str(hrmnc,'%d') 'h' num2str(nrsamplespulse) 'p.bin'];
end
savetable(Y,fname);
%% Miami Vice
% Add parameters for NIGHTTRAIN and save
MTXvel              = -[0 diff(Y(5:end)./Amp2Tab)]*SampleFreq;
saverotations(MTXvel,fname,nrntpulse);

%% Graphics Position
Y                   = Y(5:end);
subplot(231)
plot(X,Y);
hold on;
plot(X,Y2,'r');
xaxis([0 40]);
xlabel('Time');
ylabel('Amplitude');

tekst(0.1,0.9,['Maximum Amplitude =' num2str(max(abs(Y))/1000,4)]);
subplot(232)
window              = 2^15;
overlap             = 0;
[Pyy,f]             = pwelch(Y,window,overlap,[],SampleFreq);
Pyy                 = sqrt(Pyy*SampleFreq);
plot(f,Pyy,'-');

subplot(233)
sel = f<fcutoff;
plot(f(sel),Pyy(sel),'.-');
ylabel('Power (units)');


%% Graphics velocity
Y                   = [0 diff(Y)]*SampleFreq;
subplot(234)
plot(X,Y);
xaxis([0 Duration]);
xlabel('Time');
ylabel('Velocity');
xaxis([0 40]);

tekst(0.1,0.9,['Maximum Amplitude =' num2str(max(abs(Y))/1000,4)]);

subplot(235)
window      = 2^15;
overlap     = 0;
[Pyy,f]     = pwelch(Y,window,overlap,[],SampleFreq);
Pyy         = sqrt(Pyy*SampleFreq);
plot(f,Pyy,'-');

subplot(236)
sel = f<fcutoff;
plot(f(sel),Pyy(sel),'.-');
ylabel('Power (units)');



function Sig = envelope(Sig, NEnv)
SigLen      = length(Sig);
if (SigLen < 2*NEnv)
    disp ('-- ERROR: Envelope length greater than signal');
else
    Env1    = ( sin(0.5*pi*(0:NEnv)/NEnv) ).^2;
    Env2    = Env1((NEnv+1):-1:1);
    head    = 1:(NEnv+1);
    tail    = (SigLen-NEnv):SigLen;
    for i   = 1:size(Sig,2)
        Sig(head,i) = Env1' .* Sig(head,i);
        Sig(tail,i) = Env2' .* Sig(tail,i);
    end;
end;
