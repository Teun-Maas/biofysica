
h=m2c_new_tcp('192.168.1.10', 9094);


data=m2c_RD(h, 1, 'D', 1000, 2999);
data
m2c_close(h);
