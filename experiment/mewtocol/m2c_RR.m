function data = m2c_RR(handle, station_nr, start_addr, end_addr)

   msg = sprintf('0%.3d%.3d', start_addr, end_addr);
   %disp(msg);
   m2c_send(handle, station_nr, 'RR', msg);
   %disp('after m2c_send');
   data = m2c_recv(handle);
end
