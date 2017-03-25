function ic_predictcall(strf, calldatafile, h, kfreq, dt, df, loglin)
%
% IC_PREDICTCALL(strf, calldatafile, h, kfreq, dt, df, loglin);
%
%   strf: from strf1 or strf2
%   calldatafile: file with responses to calls
%	h: structure with timetraces of calls (load timetraces)
%   kfreq: low BF: 1, medium BF: 2, high BF: 3
%   dt: time bin (default: 12.5 ms);
%   df : frequency bin (default: 0.25 oct);
%   loglin: 0 for lin, 1 for log (default);
%

% Huib Versnel/John van Opstal
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<7
   loglin = 1; % use logaritmic scale for sound magnitude (power (dB))
end;
if nargin<6
   df = 0.25; % (oct) = 1/(density step)
end;
if nargin<5
   dt = 12.5; % (ms) = 1/(velocity step)
end;
Nf		= size(strf,1);
Nt		= size(strf,2);
tms		= (0:Nt-1)*dt;
freq	= (0:Nf-1)*df; % oct

%% Shift STRF-peak to centre of image
% this ensures that convolution with sound spectrogram will be in the best
% frequency range of the STRF
[~,k]	= pa_spk_strfshiftmax(strf);
if kfreq == 1, % QUE?
    k = min(k,0);
elseif kfreq==3,
    k = max(k,4);
end
strfsh	= pa_spk_shiftmatc(strf,k); % shifted strf

%% Graphics
ic_plotstrf(strf,strfsh,tms,freq)

%% Vocalizations
figure(102);
f1		= 250*2^((1-k)*0.25+2.5*(kfreq-1));
f2		= f1*2^2.25;
specs	= ic_specvocal(h,f1,f2,df,dt);

%% Prediction
ic_predcallresp(strfsh,calldatafile,specs,dt,df,f1,loglin);

function ic_plotstrf(strf,strfsh,tms,freq)
% IC_PLOTSTRF
%
% Plot shifted, unshifted, smoothed and unsmoothed STRFs

mxscal	= max(abs(strf(:))); % scaling factor for graphic


figure(101);
clf;
subplot(221);
imagesc(tms,freq,strf,[-mxscal mxscal]);
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',freq(1:2:end),'XTick',tms(1:2:end));
axis square;
colorbar;
title('Unshifted STRF');
ylabel('Frequency (oct)')

subplot(222)
x		= tms;
y		= freq;
z		= strf;
[~,m]	= max(max(z,[],1));
mxscal	= max(abs(z(:)));
XI		= linspace(min(x),max(x),500);
YI		= linspace(min(y),max(y),500);
ZI		= interp2(x,y',z,XI,YI','cubic');
imagesc(XI,YI,ZI);
shading flat;
colorbar;
axis square
caxis([-mxscal,mxscal]);
title('Unshifted, cubic-smoothed STRF');
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',freq(1:2:end),'XTick',tms(1:2:end));

subplot(223);
imagesc(tms, freq,strfsh,[-mxscal mxscal]);
xlabel('time (ms)');
title('Shifted STRF');
ylabel('Frequency (oct)')
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',freq(1:2:end),'XTick',tms(1:2:end));
axis square
colorbar;

subplot(224)
x		= tms;
y		= freq;
z		= strfsh;
[~,m]	= max(max(z,[],1));
mxscal	= max(abs(z(:)));
XI		= linspace(min(x),max(x),500);
YI		= linspace(min(y),max(y),500);
ZI		= interp2(x,y',z,XI,YI','cubic');
imagesc(XI,YI,ZI);
shading flat;
colorbar;
axis square
caxis([-mxscal,mxscal]);
xlabel('time (ms)');
title('Shifted, cubic-smoothed STRF');
set(gca,'YDir','normal','TickDir','out');
set(gca,'YTick',freq(1:2:end),'XTick',tms(1:2:end));


