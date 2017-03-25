function [B, k2] = pa_spk_shiftmatc(A, k1)
% B = IC_SHIFTMATC(A,K)
%
% Shift matrix peak to centre
%

% Huib Versnel/John van Opstal
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

lcol	= size(A,1);
k2		= mod(k1,lcol);
B		= [A(end-k2+1:end,:); A(1:end-k2,:)];
