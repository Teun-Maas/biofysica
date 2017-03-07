function data=m2c_toint(hexstr)

   s = reshape(hexstr, 4, []);
   s = s([3 4 1 2], :);
   s = s';
   x=hex2dec(s);
   
   data=typecast(uint16(x), 'int16');
end
