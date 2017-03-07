
h=m2c_new_tcp('192.168.1.10', 9094);


data=m2c_RR(h, 1, 0, 495);
data
m2c_close(h);
