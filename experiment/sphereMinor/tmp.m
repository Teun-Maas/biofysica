close all
clearvars


Fs = 48828.125;
nfft = 2^11

a = Fs/2/nfft


f = linspace(0,round(Fs/2),nfft);
f = ones(1,nfft);
f = cumsum(f);
f = f./nfft*round(Fs/2);
t = 1:nfft;
t = t/Fs;
y = sin(2*pi*f.*t);
% y = repmat(y,1,20);

figure(1)
plot(y)
hold on

whos y
% soundsc(y,Fs)

[f,a,p] = getpower(y,Fs);

p = unwrap(p);
p = p-p(1);
figure(2)
subplot(121)
plot(f,a)
verline(Fs/2);
hold on

subplot(122)
plot(f,p)
hold on

Magnitude = ones(nfft/2,1);
Nbin       = length(Magnitude);
TotalPower = sum(Magnitude.^2)
NormFactor = 1.0/TotalPower;
TwoPi      = 2*pi;
Phi        = 0.0;
Power      = 0.0;
Spectrum   = zeros (1, Nbin);
for j=1:Nbin
	Spectrum(j) = Magnitude(j) * exp (1i*Phi);
	Power = Power + NormFactor*Magnitude(j).^2;
	Phi   = Phi - TwoPi*Power;
	P(j) = Phi;
end



Spectrum = [Spectrum -conj([Spectrum(1) Spectrum((Nbin):-1:2)])];
Signal   = imag(ifft(Spectrum));
whos Signal
Signal = Signal/max(Signal);
figure(1)
plot(Signal)


% soundsc(Signal,Fs)

[f,a,p] = getpower(Signal,Fs);

p = unwrap(p);
p = p-p(1);
figure(2)
subplot(121)
plot(f,a)
verline(Fs/2);

subplot(122)
plot(f,p)

