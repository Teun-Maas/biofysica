GW/20140417

BOUWEN MEX_MODBUS ONDER WINDOWS XP 32BIT


========================================


gebruik de libmodbus-5.dll van libmodbus-3.0.5

Deze kun je met gnumex direct linken aan mex_modbus.cc, alleen

wil mex.m in matlab dat hij libmodbus-5.lib heet. Even kopieren of renamen dus.


symlink maken in de huidige directory naar de include files:

ln -s ../Modbus-windows/libmodbus-3.0.5/src

zodoende kan #include <modbus/modbus.h> werken

in matlab gnumex installeren

in matlab bouwen

mex -I. -L. mex_modbus.cc libmodbus-5.lib


voila!



TIMINGS GEMETEN VAN MATLAB TOT UITGANGEN PLC

============================================



Timing van de functie modbus_write_registers 

--------------------------------------------


Seriele verbinding FPWEB2 naar PLC is ingesteld op 119200 baud

 1 word  duurt 12ms
16 words duren 23ms


15 words netto in 11ms => 1 word duurt netto 0.733ms

overhead per write command is ongeveer 11.3ms




Freerun timing
--------------


in Matlab duurt een loop van leds.set() en leds.flush() gemiddeld 23ms per iteratie

Deze waarde is zowel in MacOS als in WindowsXP gemeten.



Gate timing
-----------


Tijd van gate high naar output actief is altijd < 3ms



Trigger timing
--------------


Tijd van positieve flank op de trigger ingang naar output actief is altijd < 3ms



LED PWM frequentie (helderheidsregeling)

----------------------------------------
Deze is ingesteld op 1kHz
