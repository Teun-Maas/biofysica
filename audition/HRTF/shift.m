function y = shift (x,n)

if (n == 0)
  y = x;
elseif (n < 0)
  n = abs(n);
  y = x([(1+n):end, 1:n],:);
else % n > 0
  y = x([(end-n+1):end 1:(end-n)],:);
end;
