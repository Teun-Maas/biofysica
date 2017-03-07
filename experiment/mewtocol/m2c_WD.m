function result = m2c_WD(handle, station_nr, mem_area_code, start_addr, data)

% size(data);
   end_addr = start_addr + length(data) -1;
   s = m2c_tohex(data);

   msg = sprintf('%c%.5d%.5d%s', mem_area_code(1), start_addr, end_addr, s);
%   disp(msg);
% size(s)
% size(msg)
   m2c_send(handle, station_nr, 'WD', msg);
%   disp('after m2c_send');

   result = m2c_recv(handle);
end
