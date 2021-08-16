
plc_addr='192.168.1.10';
fn='defaultprofiles.mat';

map=m2c_memory_map('PLC Globale Variabelen.xlsx');

t1a=map.Table_1A;

delete(map);

plc=m2c_plc(plc_addr);


Table_1A=plc.IEC_read(t1a,0,2000);
plc.IEC_write(t1a, Table_1A);
delete(plc);

if ~exist(fn, 'file')
    save(fn, ...
         'Table_1A', 'Table_1B', 'Table_1C', ...
         'Table_2A', 'Table_2B', 'Table_2C');
else
    disp('not saved, file exists');
end
