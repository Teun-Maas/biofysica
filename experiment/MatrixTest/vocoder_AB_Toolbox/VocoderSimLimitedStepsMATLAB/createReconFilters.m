function g = createReconFilters(reconStyle,reconSlope,nFFT,FS,iStartArray,iEndArray)
% createReconFilters        Create reconstruction filters.
%
% Usage:
%   g = createReconFilters(reconStyle,reconSlope,nFFT,FS,iStartArray,iEndArray)
%
% Arguments:
%   reconStyle  (0) trapezoidal (the only one supported currently)
%   reconSlope  Slope (dB/octave)
%   nFFT        Number of FFT points
%   FS          Sampling frequency
%   iStartArray Where each filter should start
%   iEndArray   Where each filter should end

% Compute the frequencies (in octaves)
fBins = log((FS/nFFT * (0:(nFFT/2))+1e-50))/log(2);

nChans = length(iStartArray);

g = -Inf*ones(nFFT/2+1,nChans);

for iChan = 1:nChans
    iStart = iStartArray(iChan);
    iEnd = iEndArray(iChan);
    
    switch(reconStyle)
        case 0  % Trapezoidal
            g(iStart:iEnd,iChan) = 0;
            g(2:iStart-1,iChan) = -reconSlope * (fBins(iStart) - fBins(2:iStart-1));
            g(iEnd+1:(nFFT/2-1),iChan) = reconSlope * (fBins(iEnd) - fBins(iEnd+1:(nFFT/2-1)));
            g(:,iChan) = 10.^(g(:,iChan)/20);
        case 1  % Triangular
            iMid = floor( (iStart+iEnd)/2);
            g(iMid,iChan) = 0;
            g(2:iMid-1,iChan) = -reconSlope * (fBins(iMid) - fBins(2:iMid-1));
            g(iMid+1:(nFFT/2-1),iChan) = reconSlope * (fBins(iMid) - fBins(iMid+1:(nFFT/2-1)));
            g(:,iChan) = 10.^(g(:,iChan)/20);
        case 2  % Triangular (in linear scale)
            iMid = floor( (iStart+iEnd)/2);
            g(iMid,iChan) = 0;
            g(2:iMid-1,iChan) = -10^(reconSlope/20) * (fBins(iMid) - fBins(2:iMid-1));
            g(iMid+1:(nFFT/2-1),iChan) = 10^(reconSlope/20) * (fBins(iMid) - fBins(iMid+1:(nFFT/2-1)));
            g(:,iChan) = g(:,iChan) + 1;
            g(find(g(:,iChan)<0),iChan) = 0;
        case 3  % Triangular quadratic
            iMid = floor( (iStart+iEnd)/2);
            g(iMid,iChan) = 0;
            g(2:iMid-1,iChan) = -reconSlope * (fBins(iMid) - fBins(2:iMid-1));
            g(iMid+1:(nFFT/2-1),iChan) = reconSlope * (fBins(iMid) - fBins(iMid+1:(nFFT/2-1)));
            g(:,iChan) = 10.^( - abs(g(:,iChan)/20) .^4 );
        case 4 % Square with equivalent 18 dB width
            iMid = floor( (iStart+iEnd)/2);
            g(iMid,iChan) = 0;
            g(2:iMid-1,iChan) = -reconSlope * (fBins(iMid) - fBins(2:iMid-1));
            g(iMid+1:(nFFT/2-1),iChan) = reconSlope * (fBins(iMid) - fBins(iMid+1:(nFFT/2-1)));

            g(find(g(:,iChan)<-18), iChan) = -1000;
            g(find(g(:,iChan)>=-18), iChan) = 0;
            
            g(:,iChan) = 10.^(g(:,iChan)/20);
            
        otherwise
            error(['reconStyle = ' num2str(reconStyle) ' not supported.']);
    end
end
