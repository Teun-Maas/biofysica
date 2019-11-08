function y = replacexy(x,new,old)


n = numel(new);
y = x;
for ii = 1:n
    y(x == old(ii)) = new(ii);
end
