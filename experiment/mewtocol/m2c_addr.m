function [ addr_out, bit_out ] = m2c_addr(addr_in, varargin)
%
% M2C_ADDR -- convert variable format address argument to numeric address
%
   if nargin == 1

      if isstr(addr_in)
         x = hex2dec(addr_in);
      else
         addr_out = addr_in;
      end
      bit_out = int16(mod(x,16));
      addr_out = int16((x-bit_out)./16); 

   elseif nargin == 2

      if isstr(addr_in)
         addr_out = hex2dec(addr_in);
      else
         addr_out = addr_in;
      end
      bit_in = varargin{1};
      if isstr(bit_in)
         bit_out = hex2dec(bit_in);
      else
         bit_out = bit_in;
      end

   else
      error('Expected one or two arguments');
   end

end
