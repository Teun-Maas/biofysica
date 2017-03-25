function mc=ic_meancorrelation(C);
%
%	function mc=meancorrelation(C);
%	
%	C: correlation coefficient matrix
%	the average of the coefficients is calculated
%
n=size(C,1);
s=sum(1:n-1);		% number of coefficients
sc=sum(sum(tril(C,-1)));
mc=sc/s;
