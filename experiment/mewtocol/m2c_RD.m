function data = m2c_RD(handle, station_nr, mem_area_code, start_addr, end_addr)

   msg = sprintf('%c%.5d%.5d', mem_area_code(1), start_addr, end_addr);
   %disp(msg);
   m2c_send(handle, station_nr, 'RD', msg);
   %disp('after m2c_send');
   data = m2c_recv(handle);
end
