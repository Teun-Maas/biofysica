% Calculates the cross-correlation coefficients for different time lags.
% 
% format: Y = MXCORR(X1,X2,LAGS);
% 
% X1 and X2 are vectors of equal size, containing the data. Use X2=X1 for autocorrelations. 
% LAGS is a vector containing the shifts between X1 and X2 in number of steps.
% MXCORR does not suffer from NaN's or an offset in the data (mean(X)~=0).
% 
% (c) Beerend Winkelman 2004