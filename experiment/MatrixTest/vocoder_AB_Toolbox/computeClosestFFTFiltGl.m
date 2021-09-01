function [fR,iStartArray,iEndArray] = computeClosestFFTFiltGl(fL, fH, FFTSTEP)
% computeClosestFFTFiltGl
%
% Optimize FFT fit globally by picking the ultimate 
%
nChans = length(fL);
fC = sqrt(fL.*fH);

errs = zeros(size(nChans));
for iAnchor = 1:(nChans-1)
    [fR,iStartArray,iEndArray] = computeClosestFFTFilt1(fL, fH, FFTSTEP, iAnchor);
    fR = (iStartArray+iEndArray)/2*FFTSTEP;
    errs(iAnchor)=sum(abs(log(fR)-log(fC)));
end
[ignore,iAnchor] = min(errs);
[fR,iStartArray,iEndArray] = computeClosestFFTFilt1(fL, fH, FFTSTEP, iAnchor);
disp(iAnchor);
