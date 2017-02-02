function Mag = remecho(Mag,Fs,grph)
% MAG = REMCHO(MAG,FS)
%
% Remove Echos from Magnitude response
%
% See also READHRTF, SWEEP2SPEC, FFTZERO, ALIGNIR, ALIGNIR2

%% Initialization
Flo                     = 400;		% Low Frequecny
Fhi                     = 19000;	% High Frequency
Fn                      = Fs/2;		% Nyquist Frequency
dt                      = 1/Fs;		% Time step in s
ndelay                  = 50;		% original = 50
if nargin<3
	grph				= 0;
end

%% Smooth
% Zero Magnitude response before Flo and after Fhi
Mag                     = fftzero(Mag, Flo/Fn, Fhi/Fn);
if grph
	magsmooth = hsmooth(Mag);
	plot(magsmooth);
	set(gca,'XScale','log');
end
%% Impuls respons
ir = ifft(Mag);

if grph
figure;
end
	if grph
		t = 1:512;
		t = 1000*t/Fs;
		m = ir(1:512,:);
		% 	m = 20*log10(m./mean(m));
		plot(t,m,'-','LineWidth',2);
		xlabel('Time (ms)');
		ylabel('Amplitude');
		drawnow
	end

% Setting the delay to 1 ms (coarse aligning)
ir                      = alignir(ir, ndelay);
% Fine alignment
ir                      = alignir2 (ir);
% Suppress reflections beyond 3 msec
Tsup                    = 0.003; % s
nsup                    = round(Tsup/dt);
ir(nsup:end-nsup+2,:)   = 0;
% Back to the Future/Frequency Domain
Mag                     = abs(fft(ir));