
h=m2c_new_tcp('192.168.1.10', 9094);
tic
idata=m2c_readint(h, 1, 1000, 2000);
toc
m2c_close(h);

pos=0.1*cumsum(idata);
times=0.1*[0:size(idata)-1];
plot(times,pos);

