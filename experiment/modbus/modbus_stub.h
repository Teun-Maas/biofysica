// debug helpers
typedef int modbus_t;

void modbus_close(modbus_t *)
{
}

void modbus_free(modbus_t *)
{
}

modbus_t* modbus_new_tcp(const char *ip_address, int port)
{
}

const char *modbus_strerror(int)
{
}

int modbus_connect(modbus_t *ctx)
{
}

int modbus_write_register(modbus_t *ctx, int reg_addr, int value)
{
}

int modbus_read_registers(modbus_t *ctx, int addr, int nb, uint16_t *dest)
{
}
