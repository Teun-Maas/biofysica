
h=m2c_new_tcp('192.168.1.10', 9094);

start=1000;
data=m2c_RD(h, 1, 'D', start, start+1999);

m2c_close(h);
idata=m2c_toint(data);
pos=0.1*cumsum(idata);
times=0.1*[0:size(idata)-1];
plot(times,pos);

