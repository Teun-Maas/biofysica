function Mag2 = pa_fftzero (Mag, flo,fhi)
% MAG = FFTZERO(MAG,FLO,FHI)
%
% Set Magnitude response below FLO and above FHI to zero
%
% 

NFFT    = size(Mag,1); % Number of Fast Fourier Transform samples
Nf      = round(NFFT/2); % Nyquist
Nlo     = round(flo*Nf); % Low sample
Nhi     = round(fhi*Nf); % High sample
posix   = 1 + (1:Nf); 
negix   = NFFT - (1:Nf) + 1;

dc      = Mag(1,:); % Do not remove DC
pos     = Mag(posix,:); % Magnitude from 2 to NFFT/2 
neg     = Mag(negix,:); % Magnitude from NFFT2/2+1 to NFFT
Ne      = 30;
Lo      = [Nlo (Nlo+Ne-1)];
Hi      = [(Nhi-Ne+1) Nhi];
pos     = pa_rampup(pa_rampdown(pos,Hi),Lo);
neg     = pa_rampup(pa_rampdown(neg,Hi),Lo);
Mag2    = [dc; pos; neg((end-1):-1:1,:)];

