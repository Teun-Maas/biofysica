function result = m2c_tohex(data)
% MEWTOCOL_TOHEX - convert data to byte swapped hex character string
%
% Syntax
%    hexstr = m2c_tohex(data);
%
   if isa(data,'int16')
       d16=typecast(data,'uint16');
       s=dec2hex(d16,4);
   elseif isa(data,'uint16')
       s=dec2hex(data,4);
   else
       error('m2c_tohex: unsupported data type');
   end
   s = s(:,[3 4 1 2]);
   result = reshape(s', 1, []);
end
