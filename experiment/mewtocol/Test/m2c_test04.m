
h=m2c_new_tcp('192.168.1.10', 9094);
addr=200;
for i=1:20
   PV_Table_Index=m2c_readint(h, 1, addr, 6);
   fprintf('%6d ', PV_Table_Index);
   fprintf('\n');
end;

m2c_close(h);

