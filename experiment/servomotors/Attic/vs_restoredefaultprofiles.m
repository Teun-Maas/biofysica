
plc_addr='192.168.1.10';
fn='defaultprofiles-panasonic-merfeld.mat';

map=m2c_memory_map('PLC Globale Variabelen.xlsx');

t1a=map.Table_1A;
t1b=map.Table_1B;
t1c=map.Table_1C;

t2a=map.Table_2A;
t2b=map.Table_2B;
t2c=map.Table_2C;

delete(map);

plc=m2c_plc(plc_addr);


load(fn);
plc.IEC_write(t1a,Table_1A,0,2000);
%Table_1B=plc.IEC_read(t1b,0,2000);
%Table_1C=plc.IEC_read(t1c,0,2000);
%Table_2A=plc.IEC_read(t2a,0,2000);
%Table_2B=plc.IEC_read(t2b,0,2000);
%Table_2C=plc.IEC_read(t2c,0,2000);
Table_1A_r = plc.IEC_read(t1a,0,2000);

delete(plc);

