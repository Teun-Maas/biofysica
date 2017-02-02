CFLAGS=-Wall -I/usr/local/include 
LDFLAGS=-L/usr/local/lib

LDLIBS=-lmodbus

mex_modbus: mex_modbus.cc
	/Applications/MATLAB_R2013b.app/bin/mex -v -I/usr/local/include -L /usr/local/lib -lmodbus mex_modbus.cc 

