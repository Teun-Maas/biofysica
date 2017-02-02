function y = pa_gsmooth(x,fs,sd,sdextra)
% Y = PA_GSMOOTH(X,FS,SD,SDEXTRA)
%
% Gaussian smoothing
%
% x  = data (samples)
% y  = smoothed data
% fs = sample frequency data (samples/s = Hz)
% sd = standard deviation gauss (s/SD)     = fs*sd samplepoints/SD
%
% sdextra = number of sd's added to data edges (5 = default)
%
% 2011 Marc van Wanrooij
% e-mail: m.vanwanrooij@neural-code.com

if nargin<4, sdextra=5; end;


% Switch data (why???)
x           = x(:)';
% Extra sample points?
%
nextra      = round(sdextra*sd*fs);
% Length(x)
nx          = length(x);
% Number of fft-points + 2 times extra sd's
nfft        = nx+2*nextra;

% determine frequency for each point
if rem(nfft,2),    % nfft odd
    frq     = [0:(nfft-1)/2 (nfft-1)/2:-1:1];
else
    frq     = [0:nfft/2 nfft/2-1:-1:1];
end

% Add extra sd
x           = [x(nextra+1:-1:2) x x(end:-1:end-nextra+1)]; % mirror edges to go around boundary effects

% Determine gaussian
g           = normpdfun(frq,0,sd*fs);
g           = g./sum(g);
% plot(frq,fft(g))
% hold on
% plot(frq,fft(x),'r')
% return
y           = real(ifft(fft(x).*fft(g)));
y           = y(nextra+1:end-nextra);

