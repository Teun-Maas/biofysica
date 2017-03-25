function data = pa_spk_ripple2strf(Spike,varargin)
% STRF = PA_SPK_RIPPLE2STRF(SPIKE)
%
% Determine STRF from SPIKE-structure (obtained by PA_SPK_READSRC). The
% STRF-structure contains the following fields:
% - magnitude
% - phase
% - density
% - velocity
% - norm
% - moddepth
% - shift
% - strf
% - time
% - frequency
%
% PA_SPK_RIPPLESTRF checks for the ripple parameters in SPIKE (velocity,
% density, depth) but nor for extraneous parameters such as duration of
% static part. This means that you have to divide the experiments in
% different SPIKE structures.
%
% PA_SPK_RIPPLE2STRF(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
%	'shift'	- shift the spectro-temporal receptive field by X octaves.
%	Default: 0 (no shift
%	'fs'	- sample frequency of ripple stimuli. Default: 50000 (Hz)
%	'nfft'	- number of FFT points for stimulus generation. Default: 2^14
%	'nbin'
%
% See also PA_SPK_READSRC, PA_SPK_PLOTSTRF

% 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com


%% Initialization
fs        = pa_keyval('fs',varargin);  % sample frequency of ripple stimulus
if isempty(fs)
	fs				= 50000;
end
nfft        = pa_keyval('nfft',varargin); % number of FFT points for stimulus generation
if isempty(nfft)
	nfft			= 2^14;
end
Nbin        = pa_keyval('nbin',varargin); % number of FFT points for stimulus generation
if isempty(Nbin)
	Nbin			= 2^5;
end
starttime        = pa_keyval('start',varargin); % Remove first X ms from analysis
if isempty(starttime)
	starttime		= 100; % ms
end
endtime        = pa_keyval('end',varargin); % End at X ms for analysis
if isempty(endtime)
	endtime		= 1000; % ms
end
df				= fs/nfft/2;	% frequency binning in ripple spectrum
x				= linspace(1/Nbin,1,Nbin);

hfshift        = pa_keyval('shift',varargin); % Do we need to shift the STRF by 1 octave
if isempty(hfshift)
	hfshift			= 0;
end
ncomp        = pa_keyval('comp',varargin); % Do we need to shift the STRF by 1 octave
if isempty(ncomp)
	ncomp			= 1;
end


%% Obtain firing rate during different periods
A	= [Spike.spiketime];
A	= A-300; % Remove the default 300 ms silent period
t	= [Spike.trial];
Dstat		= NaN(size(Spike));
Drand		= Dstat;
Velocity	= NaN(size(Spike));
Density		= Dstat;
MD			= Dstat;
for ii = 1:length(Spike)
	Dstat(ii)		= Spike(ii).stimvalues(3);
	Drand(ii)		= Spike(ii).stimvalues(4);
	Velocity(ii)	= Spike(ii).stimvalues(5);
	Density(ii)		= Spike(ii).stimvalues(6);
	MD(ii)			= Spike(ii).stimvalues(7);
end
Dstat		= Dstat(t);
A			= A-Dstat;
% Some rounding, or else Matlab will not recognize certain values correctly
% due to rounding precision
Velocity	= round(Velocity*1000)/1000;
Density		= round(Density*1000)/1000;
MD			= round(MD);
vel			= Velocity;
dens		= Density;
moddepth	= MD;
% Set Velocity, Density and Depth in correct order
Velocity	= Velocity(t);
Density		= Density(t);
MD			= MD(t);
% And get unique values
uV			= unique(Velocity);
uD			= unique(Density);
uM			= unique(MD);

%% Get STRF
M	= squeeze(NaN([numel(uV) numel(uD) numel(uM)]));
P	= M;
Nrm = M;
for ii = 1:numel(uV)
	for jj = 1:numel(uD)
		for kk = 1:numel(uM)
			selv	= Velocity==uV(ii);
			seld	= Density==uD(jj);
			selm	= MD==uM(kk);
			sela	= A>=0 & A<unique(Drand); % Select for velocity and density
			sel		= selv & seld & selm & sela;
			Nrep	= sum(vel==uV(ii) & dens==uD(jj) & moddepth == uM(kk));
			if sum(sel)
				%% Bin
				SpikeTimes	= A(sel);
				w			= round(uV(ii)/df)*df;
				period		= 1000/w;
				Nspikes		= getbins(SpikeTimes,Nbin,x,period,starttime,endtime,Nrep);
				
				%% Fourier analysis - Take first component
				X			= fft(Nspikes,Nbin)/Nbin;
				M(ii,jj,kk)		= abs(X(ncomp+1));
				P(ii,jj,kk)		= angle(X(ncomp+1));
				binnum		= Nbin;
				nFFT		= 1+binnum/2;
				Nrm(ii,jj,kk)	= abs(X(ncomp+1))/sqrt(sum(abs(X(2:nFFT)).^2));
			end
		end
	end
end
data	= struct([]);
for ii = 1:numel(uM)
	Mag		= squeeze(M(:,:,ii));
	Phase	=  squeeze(P(:,:,ii));
	Norm	= squeeze(Nrm(:,:,ii));
	
	NO		= numel(uD);
	NW		= numel(uV);
	Whi		= max(uV);
	Wlo		= min(uV);
	Wstep	= (Whi-Wlo)/(NW-1);
	Ohi		= max(uD);
	Olo		= min(uD);
	Ostep	= (Ohi-Olo)/(NO-1);
	tms		= (0:NW*2-1)/(NW*2)/Wstep*1000;  % 10 points in strf (changed 1-'01)
	xoct	= (1:NO-1)/(NO-1)/Ostep;
	
	if ~any(isnan(Mag))
		%% STRF
		C		= Mag.*exp(1i*Phase);
		strf	= compstrf(C);
	else
		% Find fixed density and velocity
		indxFixD = find(sum(isnan(Mag))==0);
		indxFixV = find(sum(isnan(Mag'))==0);
		data(ii).fixeddensity	= uD(indxFixD);
		data(ii).fixedvelocity	= uV(indxFixV);
		
		% first determine parameters for fixed density O Omega
		Mfd		= Mag(:,indxFixD);
		Pfd		= Phase(:,indxFixD);
		
		% then determine parameters for fixed velocity w omega
		Mfv		= Mag(indxFixV,:);
		Pfv		= Phase(indxFixV,:);
		% 		% index for two crosspoints
		Pcross	= Phase(indxFixV,indxFixD); % phases at cross points
		Across	= Mag(indxFixD,indxFixD); % amplitudes at cross points
		Mt		= Mfd*Mfv;
		Pt		= Pfd*Pfv;
		
		data(ii).magextra	= Mt;
		data(ii).phasextra	= Pt;
		
		Cplx	= Mt.*exp(1i*Pt);
		CplxQ12	= Cplx*exp(-1i*Pcross)/Across;
		
		% total complex transfer function and strf
		strf = compstrf(CplxQ12);
	end
	if hfshift~=0
% 		strf = pa_spk_strfshiftmax(strf);
		kshift	= hfshift/(xoct(2)-xoct(1));
		m		= size(strf,1);
		indx	= [kshift+1:m 1:kshift];
		strf	= strf(indx,:);
	end
	%% Save data
	data(ii).strf		= strf;
	data(ii).frequency	= xoct;
	data(ii).time		= tms;
	data(ii).magnitude	= Mag;
	data(ii).phase		= Phase;
	data(ii).norm		= Norm;
	data(ii).moddepth	= uM(ii);
	data(ii).density	= uD;
	data(ii).velocity	= uV;
	data(ii).shift		= hfshift;
end

function Nspikes = getbins(SpikeTimes,Nbin,x,period,starttime,endtime,Nrep)
% GETBINS
% Bin the spike responses

%% Remove Onset Transient
% Remove all periods that have spikes in the first 100 ms
N				= ceil(starttime/period); % Number of periods
% The next part could be made better
% sel				= SpikeTimes>starttime & SpikeTimes<=N*period;
% SPKAfterOnset		= SpikeTimes(sel);

sel				= SpikeTimes>N*period;
SpikeTimes		= SpikeTimes(sel); % remove all spikes that fall within the "onset transient"

%% Use only periods that fall completely within 1000 ms
N			= floor(endtime/period);
NbinafterNperiods = round(Nbin*(endtime/period-N));
sel			= SpikeTimes>N*period;
SPKAfterNperiods	= SpikeTimes(sel);
sel			= SpikeTimes<=N*period;
SpikeTimes	= SpikeTimes(sel);

%% Make period histograms
% Express Spike Times as lying between 0 and 1 period
SpikeTimes	= mod(SpikeTimes,period);
SpikeTimes	= SpikeTimes/period;
Nspikes		= hist(SpikeTimes,x); % Bin spiketimes in Nbin bins
Nspikes		= Nbin*1000*Nspikes/Nrep/period; % Convert to firing rate


if ~isempty(SPKAfterNperiods)
	SPKAfterNperiods				= mod(SPKAfterNperiods,period);
	SPKAfterNperiods				= SPKAfterNperiods/period;
	NspikesAfterNperiods	= hist(SPKAfterNperiods,x); % Bin spiketimes in Nbin bins
	NspikesAfterNperiods	= Nbin*1000*NspikesAfterNperiods/Nrep/period; % Convert to firing rate
	Nspikes					= Nspikes+NspikesAfterNperiods;
end

% 		if ~isempty(SPKAfterOnset)
% 			SPKAfterOnset				= mod(SPKAfterOnset,period);
% 			SPKAfterOnset				= SPKAfterOnset/period;
% 			NspikesAfterOnset	= hist(SPKAfterOnset,x); % Bin spiketimes in Nbin bins
% 			NspikesAfterOnset	= Nbin*1000*NspikesAfterOnset/Nrep/period; % Convert to firing rate
% 			Nspikes					= Nspikes+NspikesAfterOnset;
% 		end

%% Normalize for N periods
Nperiods = repmat(N,1,Nbin); % N full periods
Nperiods(1:NbinafterNperiods) = Nperiods(1:NbinafterNperiods)+1; % Normalize for partial period after N full periods
Nspikes = Nspikes./Nperiods;

function strf = compstrf(cplxq12,nw,no)
%  STRF = COMPSTRF(CPLXQ12,NW,NO)
%
% Obtain spectrotemporal receptive field STRF from the complex transfer
% function for quadrants 1 and 2 CPLXQ12, for w>0, entire omega range.

% Huib Versnel, version 26-7-'00
if nargin < 3
	no=1;
end;
if nargin < 2
	nw=1;
end;
s				= size(cplxq12,2);
stat			= zeros(1,s);
cplxq34			= conj(rot90(cplxq12,2));
cplxtotal		= [cplxq34;stat;cplxq12];
cplxtot1		= conj(rot90(fftshift(cplxtotal(2:end,2:end)),2));
cplxstrf		= ifft2(cplxtot1*(nw*no));	%inverse Fourier transform
strf			= real(fliplr(cplxstrf).');

