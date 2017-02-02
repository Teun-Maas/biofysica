mb = modbus_new_tcp('131.174.140.24', 502);
modbus_connect(mb);
r=modbus_read_registers(mb, 0, 5);
modbus_close(mb);
modbus_free(mb);
