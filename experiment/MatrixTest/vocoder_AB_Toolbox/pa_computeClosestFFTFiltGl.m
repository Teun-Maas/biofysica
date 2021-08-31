function [fR,indxStart,indxEnd] = pa_computeClosestFFTFiltGl(fL, fH, FFTSTEP)
% [FR,ISTART,IEND] = PA_COMPUTECLOSESTFFTFILTGL
%
% Optimize FFT fit globally by picking the ultimate 
%

% (c) 2007 Litvak, Advanced Bionics
% 2012 Modified  by Marc van Wanrooij
nChans	= length(fL);
fC		= sqrt(fL.*fH);

errs	= NaN(size(nChans));
for iAnchor = 1:(nChans-1)
    [~,indxStart,indxEnd]	= computeClosestFFTFilt1(fL,fH,FFTSTEP,iAnchor);
    fR						= (indxStart+indxEnd)/2*FFTSTEP;
    errs(iAnchor)			= sum(abs(log(fR)-log(fC)));
end
[~,iAnchor]				= min(errs);
[fR,indxStart,indxEnd]	= computeClosestFFTFilt1(fL, fH, FFTSTEP, iAnchor);

function [fR,iStartArray,iEndArray] = computeClosestFFTFilt1(fL, fH, FFTSTEP,iAnchor)

% Parameters
if nargin < 4
    iAnchor = 5;
end

% Computation
fC					= sqrt(fL.*fH);
nChans				= length(fC);
iMidArray			= zeros(1,nChans);

% Fit on the middle
iMid				= round( fC(iAnchor)/FFTSTEP );
iMidArray(iAnchor)	= iMid;


iMid				= iMid-1;
for iChan = (iAnchor-1):-1:1
    fTarget			= log(fC(iChan));
    fCur			= log(iMid*FFTSTEP);
    fNext			= log((iMid-1)*FFTSTEP);
    
    while( abs(fTarget-fCur) > abs(fTarget-fNext) )
        iMid		= iMid - 1;
        fCur		= log(iMid*FFTSTEP);
        fNext		= log((iMid-1)*FFTSTEP);
	end
    iMidArray(iChan) = iMid;
    iMid			= iMid-1;
end

iMid	= round( fC(iAnchor)/FFTSTEP )+1;
for iChan	= (iAnchor+1):1:nChans
    fTarget		= log(fC(iChan));
    fCur		= log(iMid*FFTSTEP);
    fNext		= log((iMid+1)*FFTSTEP);
    
    while( abs(fTarget-fCur) > abs(fTarget-fNext) )
        iMid	= iMid + 1;
        fCur	= log(iMid*FFTSTEP);
        fNext	= log((iMid+1)*FFTSTEP);
    end
    
    iMidArray(iChan) = iMid;
    iMid		= iMid+1;
end

fR				= iMidArray * FFTSTEP;
iStartArray		= zeros(1,nChans);
iEndArray		= zeros(1,nChans);

% Below the anchor, bias to higher frequency; above the anchor: bias down
% Map lowering down
%         1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
nMPMax			= 100;
lowMP			= reshape([(0:nMPMax)' (0:nMPMax)']',1,(nMPMax+1)*2);
highMP			= reshape([(0:nMPMax)' (1:nMPMax+1)']',1,(nMPMax+1)*2);

for iChan = 2:iAnchor
    iStartArray(iChan) = iMidArray(iChan)-lowMP(iMidArray(iChan) - iMidArray(iChan-1));
    iEndArray(iChan)   = iMidArray(iChan)+highMP(iMidArray(iChan+1) - iMidArray(iChan));
end

for iChan = (iAnchor+1):(nChans-1)
    iStartArray(iChan) = iMidArray(iChan)-highMP(iMidArray(iChan) - iMidArray(iChan-1));
    iEndArray(iChan)   = iMidArray(iChan)+lowMP(iMidArray(iChan+1) - iMidArray(iChan));
end

iStartArray(iAnchor+1) = iEndArray(iAnchor)+1;

% First filter (assume size 1)
iStartArray(1)		= iStartArray(2)-1;
iEndArray(1)		= iStartArray(1);

% Last filter
iStartArray(nChans) = iEndArray(nChans-1)+1;
iEndArray(nChans)   = round(fH(nChans)/FFTSTEP);

% Check continuity
if(sum(iStartArray(2:nChans)-iEndArray(1:(nChans-1))-1)>0.9)
    error('Filters are not continuous');
end
