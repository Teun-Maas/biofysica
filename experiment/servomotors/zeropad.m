function y = zeropad(x, d)
% ZEROPAD(X,D) -- pad vector X with zeros so its length becomes D

%GW 20211115
assert(isvector(x),'Input argument x must be a vector');  

[m,n] = size(x);
if m==1
   nzeros = d-n;
   y=[x zeros(1,nzeros)];
else
   mzeros = d-m;
   y=[x; zeros(mzeros,1)];
end

end