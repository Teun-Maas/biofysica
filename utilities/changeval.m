function B = changeval(A, new, old)
%CHANGEVAL  Substitute values in numeric matrix
%


B = A;
n = numel(new);
for ii = 1:n
	sel		= A == old(ii);
    B(sel)	= new(ii);
end
