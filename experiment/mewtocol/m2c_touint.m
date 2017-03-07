function data=m2c_touint(hexstr)
   s = reshape(hexstr, 4, []);
   s = s([3 4 1 2], :);
   s=s';
   x=hex2dec(s);
   
   data=uint16(x);
end
