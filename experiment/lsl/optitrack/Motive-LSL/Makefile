# Makefile 
# Use this with nmake from the Developer Command Prompt for VS 2017
# 
LIBROOT=\Users\gunter\Documents\Motive-LSL-dev
NATNET_DIR=$(LIBROOT)\NatNetSDK
LSL_DIR=$(LIBROOT)\liblsl
INC=/I$(NATNET_DIR)\include /I$(LSL_DIR)\include
LIBS=liblsl32.lib NatNetLib.lib Ws2_32.lib
LINKOPTS=/link /LIBPATH:$(NATNET_DIR)\lib\x86 /LIBPATH:$(LSL_DIR)\bin $(LIBS)

motive-lsl.exe:
	cl /Ox /EHsc $(INC) motive-lsl.cpp natnetlsl.cpp $(LINKOPTS) 

clean:
	del *.obj *.exe
