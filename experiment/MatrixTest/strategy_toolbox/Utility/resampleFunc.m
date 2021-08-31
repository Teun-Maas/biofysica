function y = resampleFunc(x, par)
% INPUT:
%   - x : vector/matrix to resample
%
% FIELDS FOR PAR
%   - resRatio : resampling ratio
%   - startIdx : delay of occurence of resampled vector/matrix
%  
% OUTPUT:
%   - y : resampled vector/matrix
%
% NOTES:
%
% This script peforms a repetition of columns of a row-vector/matrix and
% inserts the result into a new vector/matrix. No interpolation is
% performed.

requiredFields = {'resRatio', 'startIdx'};
checkParamFields(par, requiredFields);

[ncols, nSamples] = size(x);
nSamplesRes = round(par.resRatio*nSamples);
idxSrc = round(((0:nSamplesRes - 1)/par.resRatio) + 0.5); 
idxTgt = (1:length(idxSrc)) + par.startIdx - 1; 
y = zeros(ncols, length(idxTgt) + par.startIdx);
y(:, idxTgt) = x(:, idxSrc);