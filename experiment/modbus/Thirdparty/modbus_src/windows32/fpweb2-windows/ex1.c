#include <modbus/modbus.h>
#include <stdio.h>

int main()
{
  int r;
  modbus_t *mb;
  uint16_t tab_reg[32];

  mb = modbus_new_tcp("131.174.140.24", 502);
  if (mb<0) {
    perror("modbus_new_tcp");
    return 1;
  }

  r=modbus_connect(mb);
  if (r<0) {
    perror("modbus_connect");
    return 1;
  }

  /* Read 5 registers from the address 0 */
  r=modbus_read_registers(mb, 0, 5, tab_reg);
  if (r<0) {
    perror("modbus_read_registers");
    return 1;
  }

  modbus_close(mb);
  modbus_free(mb);
}
