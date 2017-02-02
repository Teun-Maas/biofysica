function data = modbus_read_registers(handle, address, ndata)

data = mex_modbus('read_registers', handle, address, ndata);

end
