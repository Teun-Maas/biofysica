
h=m2c_new_tcp('192.168.1.10', 9094);
   vars=m2c_readint(h, 1, 200, 50);
   fprintf('%d \n', vars);
   fprintf('\n');

m2c_close(h);

