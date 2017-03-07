function result = m2c_WC(handle, station_nr, mem_area_code, addr, bit, value)

   if (value)
       b='1';
   else
       b='0';
   end

   % 1-bit only version
   %msg = sprintf('S%c%.3d%.1X%c', mem_area_code(1), addr, bit, b);

   nbits = 1;
   msg = sprintf('P%.1d%c%.3d%.1X%c', nbits, mem_area_code(1), addr, bit, b);
   %disp(msg);
   m2c_send(handle, station_nr, 'WC', msg);

   result = m2c_recv(handle);
end
