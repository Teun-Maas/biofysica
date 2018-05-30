function sel = iseven(X)
% SEL = PA_ISEVEN(X) 
%
% Return selection vector for parity of X,
% 0 is odd, 1 even.
%
% See also ISODD



% sel = (X - 2*floor(X/2)) == 0;
sel = ~mod(X,2);
 
