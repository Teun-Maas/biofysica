#include <modbus/modbus.h>
#include <stdio.h>
#include <unistd.h>

const int offset = 55000;

int main()
{
  int i, r;
  int value;
  modbus_t *mb;

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

  value = 0x55;
  for (i=0; i<1000; i++) {
    r=modbus_write_register(mb, offset, value);
    if (r<0) {
      perror("modbus_write_register");
      return 1;
    }
    usleep(500000);
    value = ~value;
  }

  modbus_close(mb);
  modbus_free(mb);
}
