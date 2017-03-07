
h=m2c_new_tcp('192.168.1.10', 9094);


data=m2c_RC(h, 1,  'X', '48');
data
m2c_close(h);
