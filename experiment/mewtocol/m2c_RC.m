function data = m2c_RC(handle, station_nr, mem_area_code, address, varargin)

   [addr,bit] = m2c_addr(address, varargin{:});
nbits=1;
   msg = sprintf('P%.1d%c%.3d%.1X', nbits, mem_area_code(1), addr, bit);
   %disp(msg);
   m2c_send(handle, station_nr, 'RC', msg);
   %disp('after m2c_send');
   data = m2c_recv(handle);
end
