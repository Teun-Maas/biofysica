function data = m2c_readint(handle, station_nr, start_addr, nr_regs)

   hex_data = m2c_RD(handle, station_nr, 'D', start_addr, start_addr+nr_regs-1);
   data = m2c_toint(hex_data);
end



