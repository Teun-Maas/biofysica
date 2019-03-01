close all
clearvars


%% The sweep
NFFT		= 1024;
Fstart		= 0;
Fs			= 48828.125;
Fstop		= Fs;
Nsweep		= 20;
F			= (0:(NFFT-1))*Fs/NFFT;

snd			= gensweep(NFFT,Fstart,Fstop,Fs,Nsweep);
dat			= snd;

%% throw first and last sweeps, that have ramping
Nsweep		= Nsweep-2; 
indx        = NFFT + (1:(NFFT*Nsweep)); % sample index: start at sound, +1 sweep, until Nsweep
dat			= dat(indx); 
dat			= reshape(dat, NFFT, Nsweep); % reshaping the data so that it has (usually) 1024 samples and (about) 18 sweeps

mu			= mean(dat);
dat			= dat-mu; % remove DC


nfft        = 2^(nextpow2(length(dat))); % should be 1024
Spec		= fft(dat,nfft);
M			= abs(Spec);
P			= angle(Spec);

M			= 2*M(1:nfft/2,:); % throw away half
P			= unwrap(P(1:nfft/2,:));
F			= F(1:nfft/2);


whos Spec
% 	for jj		= 1:size(data,2)
% 		Spec(:,jj)  = abs(fft(x(:,jj),nfft));
% 	end
%%
figure(1)
clf
subplot(221)
plot(snd)
nicegraph

subplot(222)
plot(dat)
nicegraph


subplot(223)
plot(F,M)
nicegraph

subplot(224)
plot(F,P)
nicegraph

