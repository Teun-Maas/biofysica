//#pragma once
//#include "natnetlsl.h"
#include <NatNetTypes.h>
#include <NatNetClient.h>



void natnetlsl_init();
void natnetlsl_write(sFrameOfMocapData* data, NatNetClient* pClient);
void natnetlsl_done();