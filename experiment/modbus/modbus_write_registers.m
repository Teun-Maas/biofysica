function modbus_write_registers(handle, address, data)

if 0
   disp(['modbus_write_registers ' num2str(address) ' (0x' dec2hex(address) ')']);
   disp(dec2bin(data,16));
end
mex_modbus('write_registers', handle, address, data);

end
